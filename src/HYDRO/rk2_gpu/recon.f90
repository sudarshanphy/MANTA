module recon_module
#include "param.h"
    implicit none
contains
  function minmod(a,b) result(c)
  !$OMP DECLARE TARGET
    real :: a, b
    real :: c, fact

    fact = sign(1.0, a) + sign(1.0, b)
    c = 0.500 * fact * min(abs(a), abs(b))

  end function minmod

  subroutine get_facevalue_tvd(solnVar, recon_plus, recon_minus, check)
        use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                            iGlo, iGhi, jGlo, jGhi, kGlo, kGhi
        implicit none
        real, dimension(1:, iGlo:, jGlo:, kGlo:), intent(in) :: solnVar
        real, dimension(NDIM, NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi), intent(inout) :: recon_minus, recon_plus
        logical, intent(in), optional :: check
        real :: fplus12, fminus12
        integer :: i, j, k, n, im, ip, jm, jp, km, kp, dir
        integer :: jpad, kpad
        real :: delp, delm, delc, dxover2, qm1, q, qp1, fact
        logical :: docheck

        real, parameter :: theta = 1.40

        docheck = .false.
        if (present(check)) then
           docheck = .true.
        end if

        ! inactive axes have no reconstruction halo
        jpad = min(1,NDIM-1)
        kpad = max(0,NDIM-2)

        do dir = 1, NDIM
        select case (dir)
           case (IAXIS)
              im = 1; ip = 1; jm = 0; jp = 0; km = 0; kp = 0
           case (JAXIS)
              im = 0; ip = 0; jm = 1; jp = 1; km = 0; kp = 0
           case (KAXIS)
              im = 0; ip = 0; jm = 0; jp = 0; km = 1; kp = 1
        end select

        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP MAP(TO: solnVar) MAP(TOFROM: recon_plus, recon_minus) &
        !$OMP FIRSTPRIVATE(im, ip, jm, jp, km, kp, dir, docheck) &
        !$OMP FIRSTPRIVATE(ilo, ihi, jlo, jhi, klo, khi, jpad, kpad) &
        !$OMP PRIVATE(qm1, q, qp1) &
        !$OMP PRIVATE(fplus12, fminus12) &
        !$OMP PRIVATE(delp, delm, delc, fact, dxover2)
        do k = klo-kpad, khi+kpad
          do j = jlo-jpad, jhi+jpad
            do i = ilo-1, ihi+1
              do n = 1, NCONSVAR_NUMBER
                  qm1 = solnVar(n,i-im, j-jm, k-km)
                  q = solnVar(n,i,j,k)
                  qp1 = solnVar(n,i+ip, j+jp, k+kp)

                  delp = qp1 - q
                  delm = q - qm1
                  delc = qp1 - qm1

                  ! tvd reconstruction: generalized minmod limiter
                  fact = min(theta*abs(delp), abs(delc)/2.0, theta*abs(delm))
                  fact = sign(fact,delc)

                  if (delp * delm < 0.0) fact = 0.0

                  dxover2 = 0.5 * fact
                  fplus12 = q + dxover2
                  fminus12 = q - dxover2
                  recon_plus(dir,n,i,j,k) = fplus12
                  recon_minus(dir,n,i,j,k) = fminus12
              end do
            end do
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        end do !dir

  end subroutine

  subroutine get_facevalue_weno(solnVar, recon_plus, recon_minus, check)
        use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                            iGlo, iGhi, jGlo, jGhi, kGlo, kGhi
        implicit none
        real, dimension(1:, iGlo:, jGlo:, kGlo:), intent(in) :: solnVar
        real, dimension(NDIM, NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi), intent(inout) :: recon_minus, recon_plus
        logical, intent(in), optional :: check
        real :: fplus12, fminus12

        real, parameter :: eps = 1.0e-6
        real, dimension(2) :: beta, alpha, alphabar, &
                                 omega, omegabar, &
                                 qpluspoly, qminuspoly
        real :: qplus12, qminus12
        integer :: i, j, k, n, im, ip, jm, jp, km, kp, dir
        integer :: jpad, kpad
        real :: delp, delm
        real :: qm1, q, qp1
        logical :: docheck
        ! coefficients for the polynomial
        real, parameter, dimension(2) :: c0 = (/1.00/2.00, 1.00/2.00/)
        real, parameter, dimension(2) :: c1 = (/-1.00/2.00, 3.00/2.00/)
        real, parameter, dimension(2) :: c_1 = (/3.00/2.00, -1.00/2.00/)
        ! linear weights
        real, parameter, dimension(2) :: dcoeff = (/2.00/3.00, 1.00/3.00/)
        real, parameter, dimension(2) :: dbarcoeff = (/1.00/3.00, 2.00/3.00/)

        docheck = .false.
        if (present(check)) then
           docheck = .true.
        end if

        ! inactive axes have no reconstruction halo
        jpad = min(1,NDIM-1)
        kpad = max(0,NDIM-2)

        do dir = 1, NDIM
        select case (dir)
           case (IAXIS)
              im = 1; ip = 1; jm = 0; jp = 0; km = 0; kp = 0
           case (JAXIS)
              im = 0; ip = 0; jm = 1; jp = 1; km = 0; kp = 0
           case (KAXIS)
              im = 0; ip = 0; jm = 0; jp = 0; km = 1; kp = 1
        end select

        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP MAP(TO: solnVar) MAP(TOFROM: recon_plus, recon_minus) &
        !$OMP FIRSTPRIVATE(im, ip, jm, jp, km, kp, dir, docheck) &
        !$OMP FIRSTPRIVATE(ilo, ihi, jlo, jhi, klo, khi, jpad, kpad) &
        !$OMP PRIVATE(qm1, q, qp1, beta) &
        !$OMP PRIVATE(alpha, alphabar) &
        !$OMP PRIVATE(omega, omegabar, qpluspoly, qminuspoly) &
        !$OMP PRIVATE(qplus12, qminus12, fplus12, fminus12) &
        !$OMP PRIVATE(delp, delm)
        do k = klo-kpad, khi+kpad
          do j = jlo-jpad, jhi+jpad
            do i = ilo-1, ihi+1
              do n = 1, NCONSVAR_NUMBER
                  qm1 = solnVar(n,i-im, j-jm, k-km)
                  q = solnVar(n,i,j,k)
                  qp1 = solnVar(n,i+ip, j+jp, k+kp)

                  beta(1) = (qp1 - q)**2
                  beta(2) = (q - qm1)**2

                  alpha(1:2) = dcoeff(1:2) / (beta(1:2) + eps)**2
                  alphabar(1:2) = dbarcoeff(1:2) / (beta(1:2) + eps)**2

                  omega(1:2) = alpha(1:2) / sum(alpha)
                  omegabar(1:2) = alphabar(1:2) / sum(alphabar)

                  qpluspoly(1) = dot_product(c0(1:2), (/q, qp1/))
                  qpluspoly(2) = dot_product(c1(1:2), (/qm1, q/))

                  qminuspoly(1) = dot_product(c_1(1:2), (/q, qp1/))
                  qminuspoly(2) = dot_product(c0(1:2), (/qm1, q/))

                  qplus12 = dot_product(omega(1:2), qpluspoly(1:2))
                  qminus12 = dot_product(omegabar(1:2), qminuspoly(1:2))

                  fplus12 = qplus12
                  fminus12 = qminus12
                  ! monotonicity check
                  if ((qplus12 - q) * (q - qminus12) < 0.000) then
                     fplus12 = q
                     fminus12 = q
                  end if
                  if (docheck) then
                     delp = qp1 - q
                     delm = q - qm1
                     ! tvd reconstruction
                     if (fplus12 < 0.000) then
                       fplus12 = q + 0.500 * minmod(delp, delm)
                     end if
                     if (fminus12 < 0.000) then
                       fminus12 = q - 0.500 * minmod(delp, delm)
                     end if
                  end if
                  recon_plus(dir,n,i,j,k) = fplus12
                  recon_minus(dir,n,i,j,k) = fminus12
              end do
            end do
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

        end do !dir

  end subroutine get_facevalue_weno

  subroutine get_facevalue_weno5z(solnVar, recon_plus, recon_minus, check)
        use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                            iGlo, iGhi, jGlo, jGhi, kGlo, kGhi
        implicit none
        real, dimension(1:, iGlo:, jGlo:, kGlo:), intent(in) :: solnVar
        real, dimension(NDIM, NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi), intent(inout) :: recon_minus, recon_plus
        logical, intent(in), optional :: check

        real :: fplus12, fminus12
        real, parameter :: eps = 1.0e-36
        real, dimension(3) :: beta, alpha, alphabar, &
                              omega, omegabar, qpluspoly, &
                              qminuspoly
        real :: qplus12, qminus12
        integer :: i, j, k, n, imm, im, ip, ipp, &
                            jmm, jm, jp, jpp, kmm, km, kp, kpp, dir
        integer :: jpad, kpad
        real :: delp, delm, absbeta_diff
        real :: qm2, qm1, q, qp1, qp2
        logical :: docheck

        ! adapted from flashx spark
        real, parameter :: n13o12 = 13.0/12.0
        !! Set WENO5 coefficients once and for all
        !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
        real, parameter, dimension(3) :: coeff1p1 = (/ 2./6., -7./6., 11./6./)
        !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
        real, parameter, dimension(3) :: coeff1p2 = (/-1./6.,  5./6.,  2./6./)
        !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
        real, parameter, dimension(3) :: coeff1p3 = (/ 2./6.,  5./6., -1./6./)
        ! linear weights for i + 1/2
        real, parameter, dimension(3) :: dcoeff = (/0.1, 0.6, 0.3/)
        !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
        real, parameter, dimension(3) :: coeff1m1 = (/-1./6.,  5./6.,  2./6./)
        !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
        real, parameter, dimension(3) :: coeff1m2 = (/ 2./6.,  5./6., -1./6./)
        !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
        real, parameter, dimension(3) :: coeff1m3 = (/ 11./6.,-7./6.,  2./6./)
        ! linear weights for i - 1/2
        real, parameter, dimension(3) :: dbarcoeff = (/0.3, 0.6, 0.1/)

        docheck = .false.
        if (present(check)) then
           docheck = .true.
        end if

        ! inactive axes have no reconstruction halo
        jpad = min(1,NDIM-1)
        kpad = max(0,NDIM-2)

        do dir = 1, NDIM
        select case (dir)
           case (IAXIS)
              imm = 2; im = 1; ip = 1; ipp = 2
              jmm = 0; jm = 0; jp = 0; jpp = 0
              kmm = 0; km = 0; kp = 0; kpp = 0
           case (JAXIS)
              imm = 0; im = 0; ip = 0; ipp = 0
              jmm = 2; jm = 1; jp = 1; jpp = 2
              kmm = 0; km = 0; kp = 0; kpp = 0
           case (KAXIS)
              imm = 0; im = 0; ip = 0; ipp = 0
              jmm = 0; jm = 0; jp = 0; jpp = 0
              kmm = 2; km = 1; kp = 1; kpp = 2
        end select

        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP MAP(TO: solnVar) MAP(TOFROM: recon_plus, recon_minus) &
        !$OMP FIRSTPRIVATE(imm, im, ip, ipp, jmm, jm, jp, jpp) &
        !$OMP FIRSTPRIVATE(kmm, km, kp, kpp, docheck, dir) &
        !$OMP FIRSTPRIVATE(ilo, ihi, jlo, jhi, klo, khi, jpad, kpad) &
        !$OMP PRIVATE(qm2, qm1, q, qp1, qp2, beta) &
        !$OMP PRIVATE(absbeta_diff, alpha, alphabar) &
        !$OMP PRIVATE(omega, omegabar, qpluspoly, qminuspoly) &
        !$OMP PRIVATE(qplus12, qminus12, fplus12, fminus12) &
        !$OMP PRIVATE(delp, delm)
        do k = klo-kpad, khi+kpad
          do j = jlo-jpad, jhi+jpad
            do i = ilo-1, ihi+1
              do n = 1, NCONSVAR_NUMBER
                  qm2 = solnVar(n,i-imm, j-jmm, k-kmm)
                  qm1 = solnVar(n,i-im, j-jm, k-km)
                  q = solnVar(n,i,j,k)
                  qp1 = solnVar(n,i+ip, j+jp, k+kp)
                  qp2 = solnVar(n,i+ipp, j+jpp, k+kpp)

                  beta(1) = n13o12*( qm2 - 2.*qm1 + q )**2 &
                     + 0.25*( qm2 - 4.*qm1 + 3.*q )**2
                  beta(2) = n13o12*( qm1 - 2.*q + qp1 )**2 &
                     + 0.25*( qm1 - qp1 )**2
                  beta(3) = n13o12*( q - 2.*qp1 + qp2 )**2 &
                     + 0.25*( 3.*q - 4.*qp1 + qp2 )**2

                  absbeta_diff = abs(beta(3) - beta(1))

                  alpha(1:3) = dcoeff(1:3) * ( 1 + absbeta_diff/(beta(1:3) + eps))
                  alphabar(1:3) = dbarcoeff(1:3) * ( 1 + absbeta_diff/(beta(1:3) + eps))

                  omega(1:3) = alpha(1:3) / sum(alpha)
                  omegabar(1:3) = alphabar(1:3) / sum(alphabar)

                  qpluspoly(1) = coeff1p1(1)*qm2 + coeff1p1(2)*qm1 + coeff1p1(3)*q
                  !dot_product(coeff1p1(1:3), (/qm2, qm1, q/))
                  qpluspoly(2) = coeff1p2(1)*qm1 + coeff1p2(2)*q + coeff1p2(3)*qp1
                  !dot_product(coeff1p2(1:3), (/qm1, q, qp1/))
                  qpluspoly(3) = coeff1p3(1)*q + coeff1p3(2)*qp1 + coeff1p3(3)*qp2
                  !dot_product(coeff1p3(1:3), (/q, qp1, qp2/))

                  qminuspoly(1) = coeff1m1(1)*qm2 + coeff1m1(2)*qm1 + coeff1m1(3)*q
                  !dot_product(coeff1m1(1:3), (/qm2, qm1, q/))
                  qminuspoly(2) = coeff1m2(1)*qm1 + coeff1m2(2)*q + coeff1m2(3)*qp1
                  !dot_product(coeff1m2(1:3), (/qm1, q, qp1/))
                  qminuspoly(3) = coeff1m3(1)*q + coeff1m3(2)*qp1 + coeff1m3(3)*qp2
                  !dot_product(coeff1m3(1:3), (/q, qp1, qp2/))

                  qplus12 = omega(1)*qpluspoly(1) + omega(2)*qpluspoly(2) + &
                            omega(3)*qpluspoly(3)
                  !dot_product(omega(1:3), qpluspoly(1:3))
                  qminus12 = omegabar(1)*qminuspoly(1) + omegabar(2)*qminuspoly(2) + &
                             omegabar(3)*qminuspoly(3)
                  !dot_product(omegabar(1:3), qminuspoly(1:3))

                  fplus12 = qplus12
                  fminus12 = qminus12
                  ! monotonicity check
                  if ((qplus12 - q) * (q - qminus12) < 0.000) then
                     fplus12 = q
                     fminus12 = q
                  end if
                  if (docheck) then
                     delp = qp1 - q
                     delm = q - qm1
                     ! tvd reconstruction
                     if (fplus12 < 0.000) then
                       fplus12 = q + 0.500 * minmod(delp, delm)
                     end if
                     if (fminus12 < 0.000) then
                       fminus12 = q - 0.500 * minmod(delp, delm)
                     end if
                  end if

                  recon_plus(dir,n,i,j,k) = fplus12
                  recon_minus(dir,n,i,j,k) = fminus12
              end do
            end do
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

        end do  !dir
  end subroutine get_facevalue_weno5z

  ! GAMER-style spatial PPM with van Leer + generalized minmod limiting.
  real function ppm_limited_slope(qm, qc, qp) result(slope)
    implicit none
    !$OMP DECLARE TARGET
    real, intent(in) :: qm, qc, qp
    real, parameter :: theta_ppm = 1.5
    real :: dl, dr, dc, dvl

    dl = qc - qm
    dr = qp - qc
    dc = 0.5 * (dl + dr)
    slope = 0.0
    if (dl*dr > 0.0) then
       dvl = 2.0*dl*dr/(dl + dr)
       slope = sign(1.0, dc) * &
               min(theta_ppm*abs(dl), theta_ppm*abs(dr), abs(dc), abs(dvl))
    end if
  end function ppm_limited_slope

  subroutine get_facevalue_ppm(solnVar, recon_plus, recon_minus, check)
    use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                            iGlo, iGhi, jGlo, jGhi, kGlo, kGhi
    implicit none
    real, dimension(1:, iGlo:, jGlo:, kGlo:), intent(in) :: solnVar
    real, dimension(NDIM, NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), &
          intent(inout) :: recon_plus, recon_minus
    logical, intent(in), optional :: check
    integer :: i, j, k, n, dir, is, js, ks
    integer :: jpad, kpad
    real :: qm2, qm1, q, qp1, qp2, slope_m, slope_c, slope_p
    real :: fplus12, fminus12, dq, q6

    ! inactive axes have no reconstruction halo
    jpad = min(1,NDIM-1)
    kpad = max(0,NDIM-2)

    do dir = 1, NDIM
       select case (dir)
       case (IAXIS)
          is = 1; js = 0; ks = 0
       case (JAXIS)
          is = 0; js = 1; ks = 0
       case (KAXIS)
          is = 0; js = 0; ks = 1
       end select

       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
       !$OMP MAP(TO: solnVar) MAP(TOFROM: recon_plus, recon_minus) &
       !$OMP FIRSTPRIVATE(is, js, ks, dir, ilo, ihi, jlo, jhi, klo, khi, jpad, kpad) &
       !$OMP PRIVATE(qm2, qm1, q, qp1, qp2, slope_m, slope_c, slope_p) &
       !$OMP PRIVATE(fplus12, fminus12, dq, q6)
        do k = klo-kpad, khi+kpad
          do j = jlo-jpad, jhi+jpad
            do i = ilo-1, ihi+1
              do n = 1, NCONSVAR_NUMBER
                qm2 = solnVar(n,i-2*is,j-2*js,k-2*ks)
                qm1 = solnVar(n,i-is,j-js,k-ks)
                q   = solnVar(n,i,j,k)
                qp1 = solnVar(n,i+is,j+js,k+ks)
                qp2 = solnVar(n,i+2*is,j+2*js,k+2*ks)

                slope_m = ppm_limited_slope(qm2, qm1, q)
                slope_c = ppm_limited_slope(qm1, q, qp1)
                slope_p = ppm_limited_slope(q, qp1, qp2)

                fminus12 = 0.5*(qm1 + q) - (slope_c - slope_m)/6.0
                fplus12  = 0.5*(q + qp1) - (slope_p - slope_c)/6.0

                dq = fplus12 - fminus12
                q6 = 6.0*(q - 0.5*(fminus12 + fplus12))
                if ((fplus12-q)*(q-fminus12) <= 0.0) then
                   fminus12 = q
                   fplus12 = q
                else if (dq*q6 > dq*dq) then
                   fminus12 = 3.0*q - 2.0*fplus12
                else if (dq*q6 < -dq*dq) then
                   fplus12 = 3.0*q - 2.0*fminus12
                end if

                fminus12 = min(max(fminus12, min(qm1,q)), max(qm1,q))
                fplus12  = min(max(fplus12, min(q,qp1)), max(q,qp1))
                recon_minus(dir,n,i,j,k) = fminus12
                recon_plus(dir,n,i,j,k) = fplus12
              end do
            end do
          end do
        end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    end do
  end subroutine get_facevalue_ppm

  subroutine recon_getcellfaces(solnVar, recon_plus, recon_minus)

        use sim_data, only: recon_method, iGlo, iGhi, jGlo, jGhi, kGlo, kGhi
        use misc_module, only: to_upper
        implicit none
        real, dimension(1:, iGlo:, jGlo:, kGlo:), intent(in) :: solnVar
        real, dimension(NDIM,NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi), intent(out) :: recon_plus, recon_minus

        integer :: i, j, k, n

        ! Initialize all output slots on the device. Only the
        ! one-cell reconstruction halo is read by get_flux.
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP MAP(FROM: recon_plus, recon_minus)
        do k = kGlo, kGhi
          do j = jGlo, jGhi
            do i = iGlo, iGhi
              do n = 1, NCONSVAR_NUMBER
                recon_plus(:,n,i,j,k) = 0.0
                recon_minus(:,n,i,j,k) = 0.0
              end do
            end do
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

        if (to_upper(trim(recon_method)) == "WENO3") then
                call get_facevalue_weno(solnVar, recon_plus, recon_minus)
        else if (to_upper(trim(recon_method)) == "TVD") then
                call get_facevalue_tvd(solnVar, recon_plus, recon_minus)
        else if (to_upper(trim(recon_method)) == "WENO5") then
                call get_facevalue_weno5z(solnVar, recon_plus, recon_minus)
        elseif (to_upper(trim(recon_method)) == "PPM") then
           call get_facevalue_ppm(solnVar, recon_plus, recon_minus)
        else
           error stop "Unknown reconstruction method"
        end if
    end subroutine recon_getcellfaces

end module recon_module
