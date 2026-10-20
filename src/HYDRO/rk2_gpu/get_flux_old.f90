module flux_module
#include "param.h"
    implicit none

    ! HLLC solver is only for pure hydro simulation
    ! HLLE can be used for both pure Hydro and MHD simulation
    ! HLLD solver is only for pure MHD simulation
contains

    subroutine riemann_solve(recon_plus, recon_minus, xF, yF, zF)
        use sim_data, only: iGlo, iGhi, jGlo, jGhi, kGlo, kGhi, flux_solver
        use misc_module, only: to_upper

        implicit none
        real(8), dimension(1:NDIM,1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(in) :: recon_plus, recon_minus
        real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(out) :: xF, yF, zF

        if (to_upper(trim(flux_solver)) == "HLLC") then
#ifndef MHD
           call riemann_solve_hllc(recon_plus, recon_minus, xF, yF, zF)
#else
           print *, "HLLC is only for pure Hydro simulation"
           stop
#endif
        else if (to_upper(trim(flux_solver)) == "HLLE") then
           call riemann_solve_hlle(recon_plus, recon_minus, xF, yF, zF)

        else if (to_upper(trim(flux_solver)) == "HLLD") then
#ifdef MHD
           call riemann_solve_hlld(recon_plus, recon_minus, xF, yF, zF)
#else
           print *, "HLLD is only for MHD simulation"
           stop
#endif
        else
           print *, "Available solvers:"
           print *, "HLLE/HLLC -> Hydro"
           print *, "HLLE/HLLD -> MHD"
           stop
        end if

    end subroutine

#ifndef MHD
    subroutine riemann_solve_hllc(recon_plus, recon_minus, xF, yF, zF)
        use sim_data, only: gamma, iGlo, iGhi, jGlo, jGhi, kGlo, kGhi, &
                            jlo, jhi, ilo, ihi, klo, khi, &
                            smalld, smalle, smallp
        implicit none
        real(8), dimension(1:NDIM,1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(in) :: recon_plus, recon_minus
        real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(out) :: xF, yF, zF

        real(8), dimension(1:NCONSVAR_NUMBER) :: FL, FR, UL, UR, Fstar, Flux, UstarL, UstarR, FstarL, FstarR
        real(8) :: dL, udL, vdL, wdL, pL
        real(8) :: dR, udR, vdR, wdR, pR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, cL, cR, HL, HR
        real(8) :: ubar, vbar, wbar, hbar, cbar
        real(8) :: qL, qR, qbar, SL, SR, SM
        real(8) :: dstarL, dstarR, pstar
        real(8) :: dustarL, dustarR, dvstarL, dvstarR
        real(8) :: dwstarL, dwstarR, estarL, estarR
        integer, dimension(3) :: n

        integer :: i, j, k, is, js, ks, dir

        do dir = 1, NDIM
        select case (dir)
        case (IAXIS)
            is = 1; js = 0; ks = 0
            n = (/1, 0, 0/)
        case (JAXIS)
            is = 0; js = 1; ks = 0
            n = (/0, 1, 0/)
        case (KAXIS)
            is = 0; js = 0; ks = 1
            n = (/0, 0, 1/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

        ! compute only the faces needed in this direction
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP MAP(TO: recon_plus, recon_minus) MAP(TOFROM: xF, yF, zF) &
        !$OMP FIRSTPRIVATE(ilo, ihi, jlo, jhi, klo, khi) &
        !$OMP MAP(TO: n, gamma, smalld, smallp, smalle) &
        !$OMP FIRSTPRIVATE(dir, is, js, ks) &
        !$OMP PRIVATE(dL, ufL, vfL, wfL, pL, eL) &
        !$OMP PRIVATE(dR, ufR, vfR, wfR, pR, eR) &
        !$OMP PRIVATE(sdR, sdL, cL, cR, qL, qR) &
        !$OMP PRIVATE(udL, udR, vdL, vdR, wdL, wdR) &
        !$OMP PRIVATE(HL, HR, ubar, vbar, wbar, Hbar, cbar, qbar) &
        !$OMP PRIVATE(dstarL, dstarR, pstar, dustarL, dustarR, dvstarL, dvstarR) &
        !$OMP PRIVATE(dwstarL, dwstarR, estarL, estarR) &
        !$OMP PRIVATE(UL, UR, SL, SR, SM, FL, FR, Fstar, Flux) &
        !$OMP PRIVATE(UstarL, UstarR, FstarL, FstarR)
        do k = klo, khi + ks
          do j = jlo, jhi + js
            do i = ilo, ihi + is
              dL  = max(recon_plus(dir,DENS_VAR,i-is,j-js,k-ks), smalld)
              ufL = recon_plus(dir,VELX_VAR,i-is,j-js,k-ks)
              vfL = recon_plus(dir,VELY_VAR,i-is,j-js,k-ks)
              wfL = recon_plus(dir,VELZ_VAR,i-is,j-js,k-ks)
              pL  = max(recon_plus(dir,PRES_VAR,i-is,j-js,k-ks), smallp)
              !eL  = recon_plus(dir,ENER_VAR,i-is,j-js,k-ks)
              ! apply eos since we can
              eL = pL/(gamma - 1.0) + 0.5 * (dL * (ufL**2 + vfL**2 + wfL**2))
              eL = max(eL, smalle)

              dR  = max(recon_minus(dir,DENS_VAR,i,j,k), smalld)
              ufR = recon_minus(dir,VELX_VAR,i,j,k)
              vfR = recon_minus(dir,VELY_VAR,i,j,k)
              wfR = recon_minus(dir,VELZ_VAR,i,j,k)
              pR  = max(recon_minus(dir,PRES_VAR,i,j,k), smallp)
              !eR  = recon_minus(dir,ENER_VAR,i,j,k)
              ! apply eos since we can
              eR = pR/(gamma - 1.0) + 0.5 * (dR * (ufR**2 + vfR**2 + wfR**2))
              eR = max(eR, smalle)

              sdR = sqrt(dR)
              sdL = sqrt(dL)

              udL = ufL * dL
              udR = ufR * dR

              vdL = vfL * dL
              vdR = vfR * dR

              wdL = wfL * dL
              wdR = wfR * dR

              HL = (eL + pL) / dL
              HR = (eR + pR) / dR

              cL = sqrt(gamma * pL / dL)
              cR = sqrt(gamma * pR / dR)

              ubar = (ufL * sdL + ufR * sdR) / (sdL + sdR)
              vbar = (vfL * sdL + vfR * sdR) / (sdL + sdR)
              wbar = (wfL * sdL + wfR * sdR) / (sdL + sdR)

              Hbar = (HL * sdL + HR * sdR) / (sdL + sdR)

              cbar = sqrt((gamma - 1.00) * (Hbar - 0.50 * (ubar**2 + vbar**2 + wbar**2)))

              ! speed in the normal direction
              qL = recon_plus(dir,VELX_VAR+dir-1,i-is,j-js,k-ks)
              qR = recon_minus(dir,VELX_VAR+dir-1,i,j,k)

              qbar = ubar*n(1) + vbar*n(2) + wbar*n(3)

              SL = min(qL - cL, qbar - cbar)
              SR = max(qR + cR, qbar + cbar)

              SM = ((dR * qR * (SR - qR) - dL * qL * (SL - qL)) + pL - pR) &
                            / (dR * (SR - qR) - dL * (SL - qL))

              FL = (/dL * qL, &
                     dL * ufL * qL + pL * n(1), &
                     dL * vfL * qL + pL * n(2), &
                     dL * wfL * qL + pL * n(3), &
                     qL * (eL + pL)/)
              FR = (/dR * qR, &
                     dR * ufR * qR + pR * n(1), &
                     dR * vfR * qR + pR * n(2), &
                     dR * wfR * qR + pR * n(3), &
                     qR * (eR + pR)/)

              dstarL = dL * (SL - qL) / (SL - SM)
              dstarR = dR * (SR - qR) / (SR - SM)

              pstar = pL + dL * (qL - SL) * (qL - SM)

              dustarL = (dL * ufL * (SL - qL) + (pstar - pL) * n(1)) / (SL - SM)
              dustarR = (dR * ufR * (SR - qR) + (pstar - pR) * n(1)) / (SR - SM)

              dvstarL = (dL * vfL * (SL - qL) + (pstar - pL) * n(2)) / (SL - SM)
              dvstarR = (dR * vfR * (SR - qR) + (pstar - pR) * n(2)) / (SR - SM)

              dwstarL = (dL * wfL * (SL - qL) + (pstar - pL) * n(3)) / (SL - SM)
              dwstarR = (dR * wfR * (SR - qR) + (pstar - pR) * n(3)) / (SR - SM)

              estarL = ((eL * (SL - qL) - pL * qL) + pstar * SM) / (SL - SM)
              estarR = ((eR * (SR - qR) - pR * qR) + pstar * SM) / (SR - SM)

              UstarL = (/dstarL, dustarL, dvstarL, dwstarL, estarL/)
              UstarR = (/dstarR, dustarR, dvstarR, dwstarR, estarR/)

              UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL/)
              UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR/)

              FstarL = FL + SL * (UstarL - UL)
              FstarR = FR + SR * (UstarR - UR)

              if (0.000 < SL) then
                  Flux = FL
              else if (SL <= 0.000 .and. 0.000 < SM) then
                  Flux = FstarL
              else if (SM <= 0.000 .and. 0.000 <= SR) then
                  Flux = FstarR
              else
                  Flux = FR
              end if

              if (dir == IAXIS) then
                xF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              else if (dir == JAXIS) then
                yF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              else if (dir == KAXIS) then
                zF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              end if

            end do
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        end do  !dir

       return

    end subroutine

#endif

    subroutine riemann_solve_hlle(recon_plus, recon_minus, xF, yF, zF)
        use sim_data, only: gamma, iGlo, iGhi, jGlo, jGhi, kGlo, kGhi, &
                            jlo, jhi, ilo, ihi, klo, khi, &
                            smalld, smalle, smallp
#ifdef MHD
        use sim_data, only: ch
#endif
        implicit none
        real(8), dimension(1:NDIM,1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(in) :: recon_plus, recon_minus
        real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(out) :: xF, yF, zF

        real(8), dimension(1:NCONSVAR_NUMBER) :: FL, FR, UL, UR, Fstar, Flux
        real(8) :: dL,  pL
        real(8) :: dR,  pR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, HL, HR, cL, cR
#ifdef MHD
        real(8) :: bxL, byL, bzL, bpL
        real(8) :: bxR, byR, bzR, bpR
        real(8) :: B2L, B2R, vBL, vBR, qBL, qBR, cfL, cfR
#endif
        real(8) :: ubar, vbar, wbar, Hbar, cbar
        real(8) :: qL, qR, qbar, SL, SR
        integer :: i, j, k, dir
        integer, dimension(3) :: n
        ! index for different dimension
        integer :: id, imn, imt1, imt2, ibn, ibt1, ibt2, ie, ips
        integer :: is, js, ks

        do dir = 1, NDIM
        id = 1; ie = 5; ips = 9
        select case (dir)
        case (IAXIS)
            imn = 2; imt1 = 3; imt2 = 4
            ibn = 6; ibt1 = 7; ibt2 = 8
            is = 1; js = 0; ks = 0
            n = (/1, 0, 0/)
        case (JAXIS)
            imn = 3; imt1 = 2; imt2 = 4
            ibn = 7; ibt1 = 6; ibt2 = 8
            is = 0; js = 1; ks = 0
            n = (/0, 1, 0/)
        case (KAXIS)
            imn = 4; imt1 = 2; imt2 = 3
            ibn = 8; ibt1 = 6; ibt2 = 7
            is = 0; js = 0; ks = 1
            n = (/0, 0, 1/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

        ! compute only the faces needed in this direction
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP MAP(TO: recon_plus, recon_minus) MAP(TOFROM: xF, yF, zF) &
        !$OMP FIRSTPRIVATE(ilo, ihi, jlo, jhi, klo, khi) &
        !$OMP MAP(TO: n, gamma, smalld, smallp, smalle) &
#ifdef MHD
        !$OMP MAP(TO: ch) &
#endif
        !$OMP FIRSTPRIVATE(dir, id, ips, ibn, is, js, ks) &
        !$OMP PRIVATE(dL, ufL, vfL, wfL, pL, eL) &
        !$OMP PRIVATE(dR, ufR, vfR, wfR, pR, eR) &
        !$OMP PRIVATE(sdR, sdL, cL, cR, qL, qR) &
#ifndef MHD
        !$OMP PRIVATE(HL, HR, ubar, vbar, wbar, Hbar, cbar, qbar) &
#endif
#ifdef MHD
        !$OMP PRIVATE(bxL, byL, bzL, bpL) &
        !$OMP PRIVATE(bxR, byR, bzR, bpR) &
        !$OMP PRIVATE(qBL, qBR, B2L, B2R, vBL, vBR, cfL, cfR) &
#endif
        !$OMP PRIVATE(UL, UR, SL, SR, FL, FR, Fstar, Flux)
        do k = klo, khi + ks
          do j = jlo, jhi + js
            do i = ilo, ihi + is
              dL  = max(recon_plus(dir,DENS_VAR,i-is,j-js,k-ks), smalld)
              ufL = recon_plus(dir,VELX_VAR,i-is,j-js,k-ks)
              vfL = recon_plus(dir,VELY_VAR,i-is,j-js,k-ks)
              wfL = recon_plus(dir,VELZ_VAR,i-is,j-js,k-ks)
              pL  = max(recon_plus(dir,PRES_VAR,i-is,j-js,k-ks), smallp)
              !eL  = recon_plus(dir,ENER_VAR,i-is,j-js,k-ks)
              ! apply eos since we can
              eL = pL/(gamma - 1.0) + 0.5 * (dL * (ufL**2 + vfL**2 + wfL**2))
              eL = max(eL, smalle)

              dR  = max(recon_minus(dir,DENS_VAR,i,j,k), smalld)
              ufR = recon_minus(dir,VELX_VAR,i,j,k)
              vfR = recon_minus(dir,VELY_VAR,i,j,k)
              wfR = recon_minus(dir,VELZ_VAR,i,j,k)
              pR  = max(recon_minus(dir,PRES_VAR,i,j,k), smallp)
              !eR  = recon_minus(dir,ENER_VAR,i,j,k)
              ! apply eos since we can
              eR = pR/(gamma - 1.0) + 0.5 * (dR * (ufR**2 + vfR**2 + wfR**2))
              eR = max(eR, smalle)
#ifdef MHD
              bxL = recon_plus(dir,BMFX_VAR,i-is,j-js,k-ks)
              byL = recon_plus(dir,BMFY_VAR,i-is,j-js,k-ks)
              bzL = recon_plus(dir,BMFZ_VAR,i-is,j-js,k-ks)
              bpL = recon_plus(dir,BPSI_VAR,i-is,j-js,k-ks)
              eL = max(eL + 0.5*(bxL*bxL + byL*byL + bzL*bzL), smalle)

              bxR = recon_minus(dir,BMFX_VAR,i,j,k)
              byR = recon_minus(dir,BMFY_VAR,i,j,k)
              bzR = recon_minus(dir,BMFZ_VAR,i,j,k)
              bpR = recon_minus(dir,BPSI_VAR,i,j,k)
              eR = max(eR + 0.5*(bxR*bxR + byR*byR + bzR*bzR), smalle)
#endif

              sdR = sqrt(dR)
              sdL = sqrt(dL)

              cL = sqrt(gamma * pL / dL)
              cR = sqrt(gamma * pR / dR)

              ! speed in the normal direction
              qL = recon_plus(dir,VELX_VAR+dir-1,i-is,j-js,k-ks)
              qR = recon_minus(dir,VELX_VAR+dir-1,i,j,k)

#ifndef MHD
              UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL/)
              UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR/)

              HL = (eL + pL) / dL
              HR = (eR + pR) / dR

              ubar = (ufL * sdL + ufR * sdR) / (sdL + sdR)
              vbar = (vfL * sdL + vfR * sdR) / (sdL + sdR)
              wbar = (wfL * sdL + wfR * sdR) / (sdL + sdR)

              Hbar = (HL * sdL + HR * sdR) / (sdL + sdR)

              cbar = sqrt((gamma - 1.00) * (Hbar - 0.50 * (ubar**2 + vbar**2 + wbar**2)))

              qbar = sum((/ubar, vbar, wbar/) * n)

              SL = min(min(qL - cL, qbar - cbar), 0.0)
              SR = max(max(qR + cR, qbar + cbar), 0.0)

              FL = (/dL * qL, &
                     dL * ufL * qL + pL * n(1), &
                     dL * vfL * qL + pL * n(2), &
                     dL * wfL * qL + pL * n(3), &
                     qL * (eL + pL)/)
              FR = (/dR * qR, &
                     dR * ufR * qR + pR * n(1), &
                     dR * vfR * qR + pR * n(2), &
                     dR * wfR * qR + pR * n(3), &
                     qR * (eR + pR)/)

              Fstar = (SR * FL - SL * FR + (SL * SR * (UR - UL))) / (SR - SL)
#endif
#ifdef MHD
              ! magnetic field value in the normal direction
              qBL = recon_plus(dir,BMFX_VAR+dir-1,i-is,j-js,k-ks)
              qBR = recon_minus(dir,BMFX_VAR+dir-1,i,j,k)

              ! B.B quantity
              B2L = bxL*bxL + byL*byL + bzL*bzL
              B2R = bxR*bxR + byR*byR + bzR*bzR

              ! V.B quantity
              vBL = ufL*bxL + vfL*byL + wfL*bzL
              vBR = ufR*bxR + vfR*byR + wfR*bzR

              ! fast magnetosonic wave speed
              ! Miyoshi eqn 67 HLLD paper
              cfL = sqrt(((gamma * pL + B2L) + &
                sqrt((gamma * pL + B2L)**2 - gamma * pL * 4.0 * qBL**2))/(2.0 * dL))

              cfR = sqrt(((gamma * pR + B2R) + &
                sqrt((gamma * pR + B2R)**2 - gamma * pR * 4.0 * qBR**2))/(2.0 * dR))

              SL = min(qL, qR) - max(cfL, cfR)
              SR = max(qL, qR) + max(cfL, cfR)

              FL = (/dL * qL, &
                     dL * ufL * qL - qBL * bxL + (0.50 * B2L + pL) * n(1) , &
                     dL * vfL * qL - qBL * byL + (0.50 * B2L + pL) * n(2) , &
                     dL * wfL * qL - qBL * bzL + (0.50 * B2L + pL) * n(3) , &
                     qL * (eL + pL + 0.50 * B2L) - qBL * vBL              , &
                     qL * bxL - qBL * ufL                                 , &
                     qL * byL - qBL * vfL                                 , &
                     qL * bzL - qBL * wfL                                 , &
                     0.0/)
              FR = (/dR * qR, &
                     dR * ufR * qR - qBR * bxR + (0.50 * B2R + pR) * n(1) , &
                     dR * vfR * qR - qBR * byR + (0.50 * B2R + pR) * n(2) , &
                     dR * wfR * qR - qBR * bzR + (0.50 * B2R + pR) * n(3) , &
                     qR * (eR + pR + 0.50 * B2R) - qBR * vBR              , &
                     qR * bxR - qBR * ufR                                 , &
                     qR * byR - qBR * vfR                                 , &
                     qR * bzR - qBR * wfR                                 , &
                     0.0/)

              UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL, bxL, byL, bzL, bpL/)
              UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR, bxR, byR, bzR, bpR/)

              Fstar = (SR * FL - SL * FR + (SL * SR * (UR - UL))) / (SR - SL)
#endif
              if (0.000 < SL) then
                  Flux = FL
              else if (SR < 0.000) then
                  Flux = FR
              else
                  Flux = Fstar
              end if
#ifdef MHD
              Flux(BPSI_CONS) = 0.50 * ((ch * ch) * (qBL + qBR) - ch * (bpR - bpL))
              Flux(BMFX_CONS+dir-1) = 0.50 * ((bpL + bpR) - ch * (qBR - qBL))
#endif
              if (dir == IAXIS) then
                xF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              else if (dir == JAXIS) then
                yF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              else if (dir == KAXIS) then
                zF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              end if

            end do
          end do
        end do
           !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
           end do !dir

           return

    end subroutine

#ifdef MHD
    subroutine riemann_solve_hlld(recon_plus, recon_minus, xF, yF, zF)
        use sim_data, only: gamma, iGlo, iGhi, jGlo, jGhi, kGlo, kGhi, &
                            jlo, jhi, ilo, ihi, klo, khi, &
                            smalld, smalle, smallp, ch
        implicit none
        real(8), dimension(1:NDIM,1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(in) :: recon_plus, recon_minus
        real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi), intent(out) :: xF, yF, zF

        real(8), dimension(1:NCONSVAR_NUMBER) :: FL, FR, UL, UR, UL1, U2, UR1, Uhll, Flux
        real(8) :: dL, pL, bxL, byL, bzL, bpL
        real(8) :: dR, pR, bxR, byR, bzR, bpR
        real(8) :: ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, cfL, cfR, ptL, ptR, pt
        real(8) :: qL, qR, SL, SR, SM, qBL, qBR, B2L, B2R, vBL, vBR
        real(8) :: f1, f2, SL1, SR1, aVL, aVR, at1, at2, vt1L, vt2L
        real(8) :: vt1R, vt2R, dsL, dsR, sdsL, sdsR, enL, enR
        real(8) :: VLibt1, VLibt2, VRibt1, VRibt2, VLips, VRips
        real(8) :: VLimt1, VRimt1, VLimt2, VRimt2, VLibn, VRibn
        real(8), parameter :: hllg_factor = 1.001d0
        integer, dimension(3) :: n
        integer :: i, j, k, dir
        ! index for different dimension
        integer :: id, imn, imt1, imt2, ibn, ibt1, ibt2, ie, ips
        integer :: is, js, ks

        do dir = 1, NDIM
        id = 1; ie = 5; ips = 9
        select case (dir)
        case (IAXIS)
            imn = 2; imt1 = 3; imt2 = 4
            ibn = 6; ibt1 = 7; ibt2 = 8
            is = 1; js = 0; ks = 0
            n = (/1, 0, 0/)
        case (JAXIS)
            imn = 3; imt1 = 2; imt2 = 4
            ibn = 7; ibt1 = 6; ibt2 = 8
            is = 0; js = 1; ks = 0
            n = (/0, 1, 0/)
        case (KAXIS)
            imn = 4; imt1 = 2; imt2 = 3
            ibn = 8; ibt1 = 6; ibt2 = 7
            is = 0; js = 0; ks = 1
            n = (/0, 0, 1/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

        ! compute only the faces needed in this direction
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP MAP(TO: recon_plus, recon_minus) MAP(TOFROM: xF, yF, zF) &
        !$OMP FIRSTPRIVATE(ilo, ihi, jlo, jhi, klo, khi) &
        !$OMP MAP(TO: n, ch, gamma, smalld, smallp, smalle) &
        !$OMP FIRSTPRIVATE(dir, id, ips, ibn, imn, imt1, imt2) &
        !$OMP FIRSTPRIVATE(ie, ibt1, ibt2, is, js, ks) &
        !$OMP PRIVATE(VLibn, VRibn, VLips, VRips) &
        !$OMP PRIVATE(VLimt1, VRimt1, VLimt2, VRimt2) &
        !$OMP PRIVATE(VLibt1, VRibt1, VLibt2, VRibt2) &
        !$OMP PRIVATE(dL, ufL, vfL, wfL, pL, eL) &
        !$OMP PRIVATE(dR, ufR, vfR, wfR, pR, eR) &
        !$OMP PRIVATE(qL, qR, ptL, ptR, dsL, dsR, aVL, aVR) &
        !$OMP PRIVATE(at1, at2, f1, f2, SL1, SR1) &
        !$OMP PRIVATE(vt1R, vt2R, vt1L, vt2L) &
        !$OMP PRIVATE(sdsL, sdsR, enL, enR, pt) &
        !$OMP PRIVATE(bxL, byL, bzL, bpL) &
        !$OMP PRIVATE(bxR, byR, bzR, bpR) &
        !$OMP PRIVATE(qBL, qBR, B2L, B2R, vBL, vBR, cfL, cfR) &
        !$OMP PRIVATE(UL, UR, SL, SR, FL, FR, Flux) &
        !$OMP PRIVATE(SM, UL1, UR1, U2, Uhll)
        do k = klo, khi + ks
          do j = jlo, jhi + js
            do i = ilo, ihi + is
              dL  = max(recon_plus(dir,DENS_VAR,i-is,j-js,k-ks), smalld)
              ufL = recon_plus(dir,VELX_VAR,i-is,j-js,k-ks)
              vfL = recon_plus(dir,VELY_VAR,i-is,j-js,k-ks)
              wfL = recon_plus(dir,VELZ_VAR,i-is,j-js,k-ks)
              pL  = max(recon_plus(dir,PRES_VAR,i-is,j-js,k-ks), smallp)
              ! apply eos since we can
              eL = pL/(gamma - 1.0) + 0.5 * (dL * (ufL**2 + vfL**2 + wfL**2))
              bxL = recon_plus(dir,BMFX_VAR,i-is,j-js,k-ks)
              byL = recon_plus(dir,BMFY_VAR,i-is,j-js,k-ks)
              bzL = recon_plus(dir,BMFZ_VAR,i-is,j-js,k-ks)
              bpL = recon_plus(dir,BPSI_VAR,i-is,j-js,k-ks)
              eL = max(eL + 0.5*(bxL*bxL + byL*byL + bzL*bzL), smalle)

              dR  = max(recon_minus(dir,DENS_VAR,i,j,k), smalld)
              ufR = recon_minus(dir,VELX_VAR,i,j,k)
              vfR = recon_minus(dir,VELY_VAR,i,j,k)
              wfR = recon_minus(dir,VELZ_VAR,i,j,k)
              pR  = max(recon_minus(dir,PRES_VAR,i,j,k), smallp)
              ! apply eos since we can
              eR = pR/(gamma - 1.0) + 0.5 * (dR * (ufR**2 + vfR**2 + wfR**2))
              bxR = recon_minus(dir,BMFX_VAR,i,j,k)
              byR = recon_minus(dir,BMFY_VAR,i,j,k)
              bzR = recon_minus(dir,BMFZ_VAR,i,j,k)
              bpR = recon_minus(dir,BPSI_VAR,i,j,k)
              eR = max(eR + 0.5*(bxR*bxR + byR*byR + bzR*bzR), smalle)

              ! some other quantities needed later:
              VLimt1 = recon_plus(dir, imt1,i-is,j-js,k-ks)
              VRimt1 = recon_minus(dir, imt1,i,j,k)
              VLimt2 = recon_plus(dir, imt2,i-is,j-js,k-ks)
              VRimt2 = recon_minus(dir, imt2,i,j,k)
              VLibt1 = recon_plus(dir, ibt1,i-is,j-js,k-ks)
              VRibt1 = recon_minus(dir, ibt1,i,j,k)
              VLibt2 = recon_plus(dir, ibt2,i-is,j-js,k-ks)
              VRibt2 = recon_minus(dir, ibt2,i,j,k)
              VLibn  = recon_plus(dir, ibn,i-is,j-js,k-ks)
              VRibn  = recon_minus(dir, ibn,i,j,k)
              VLips  = recon_plus(dir, ips,i-is,j-js,k-ks)
              VRips  = recon_minus(dir, ips,i,j,k)

              ! speed in the normal direction
              qL = recon_plus(dir,VELX_VAR+dir-1,i-is,j-js,k-ks)
              qR = recon_minus(dir,VELX_VAR+dir-1,i,j,k)

              ! magnetic field value in the normal direction
              qBL = recon_plus(dir,BMFX_VAR+dir-1,i-is,j-js,k-ks)
              qBR = recon_minus(dir,BMFX_VAR+dir-1,i,j,k)

              ! B.B quantity
              B2L = bxL*bxL + byL*byL + bzL*bzL
              B2R = bxR*bxR + byR*byR + bzR*bzR

              ! V.B quantity
              vBL = ufL*bxL + vfL*byL + wfL*bzL
              vBR = ufR*bxR + vfR*byR + wfR*bzR

              ! fast magnetosonic wave speed
              ! Miyoshi eqn 67 HLLD paper
              cfL = sqrt(((gamma * pL + B2L) + &
                sqrt(max((gamma * pL + B2L)**2 - gamma * pL * 4.0 * qBL**2,0.0)))/(2.0 * dL))

              cfR = sqrt(((gamma * pR + B2R) + &
                sqrt(max((gamma * pR + B2R)**2 - gamma * pR * 4.0 * qBR**2,0.0)))/(2.0 * dR))

              ptL = pL + 0.5e0 * B2L
              ptR = pR + 0.5e0 * B2R

              FL = (/dL * qL, &
                     dL * ufL * qL - qBL * bxL + (0.50 * B2L + pL) * n(1) , &
                     dL * vfL * qL - qBL * byL + (0.50 * B2L + pL) * n(2) , &
                     dL * wfL * qL - qBL * bzL + (0.50 * B2L + pL) * n(3) , &
                     qL * (eL + pL + 0.50 * B2L) - qBL * vBL              , &
                     qL * bxL - qBL * ufL                                 , &
                     qL * byL - qBL * vfL                                 , &
                     qL * bzL - qBL * wfL                                 , &
                     0.0/)
              FR = (/dR * qR, &
                     dR * ufR * qR - qBR * bxR + (0.50 * B2R + pR) * n(1) , &
                     dR * vfR * qR - qBR * byR + (0.50 * B2R + pR) * n(2) , &
                     dR * wfR * qR - qBR * bzR + (0.50 * B2R + pR) * n(3) , &
                     qR * (eR + pR + 0.50 * B2R) - qBR * vBR              , &
                     qR * bxR - qBR * ufR                                 , &
                     qR * byR - qBR * vfR                                 , &
                     qR * bzR - qBR * wfR                                 , &
                     0.0/)

              SL = min(qL, qR) - max(cfL, cfR)
              SR = max(qL, qR) + max(cfL, cfR)

              UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL, bxL, byL, bzL, bpL/)
              UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR, bxR, byR, bzR, bpR/)

              if (SL >= 0.000) then
                  Flux(:) = FL(:)
              else if (SR <= 0.000) then
                  Flux(:) = FR(:)
              else
                  ! Miyoshi 2005 eq 38
                  ! SM speed from HLLC solver
                  ! SM = ((dR * qR * (SR - qR) - dL * qL * (SL - qL)) + pL - pR) &
                  !                          / (dR * (SR - qR) - dL * (SL - qL))

                  ! taken from opnmhd code 2D_basic version

                  Uhll(1:8) = (SR * UR(1:8) - SL * UL(1:8) - FR(1:8) + FL(1:8)) &
                                              / (SR - SL)
                  ! entropy wave
                  SM = Uhll(imn) / Uhll(id)

                  ! total pressure: M05 eqn 23
                  pt = ptL + dL * (SL - qL) * (SM - qL)
                  ! density in star region: M05 eqn 43
                  dsL = dL * (SL - qL) / (SL - SM) !d(L*)
                  dsR = dR * (SR - qR) / (SR - SM) !d(R*)

                  ! alfven wave
                  aVL = abs(Uhll(ibn))/sqrt(dsL)
                  aVR = abs(Uhll(ibn))/sqrt(dsR)

                  ! ========== revert to HLLC-G ==========
                  !When SR*(SL*) --> SR(SL), we switch to a HLLC solver.
                  !In these limits, the HLLD denominator  rho_L (SL-uL)(SL-SM) - B_x^2,
                  !which can be written as  rho_L* (SL+SL*-2SM) (SL-SL*), becomes zero.
                  !Since the two HLLC states L* and R* are relevant to L** and R** in HLLD,
                  !we should employ HLLC-G method for logical consistency:
                  !i.e.:  vy_L* = vy_R* (HLLC-G)  <==>  vy_L** = vy_R** (HLLD).
                  !if ( .true. ) then

                  if((SL >= (SM - hllg_factor * aVL)) .or. &
                    ((SM + hllg_factor * aVR) >= SR)) then
                  f1 = 1.0e0 / (SR - SL)
                  f2 = SR * SL
                  Flux(ibt1) = f1 * (SR * FL(ibt1) - SL * FR(ibt1) + f2 * (VRibt1 - VLibt1))
                  Flux(ibt2) = f1 * (SR * FL(ibt2) - SL * FR(ibt2) + f2 * (VRibt2 - VLibt2))

                  at1 = Uhll(imt1)/Uhll(id)
                  at2 = Uhll(imt2)/Uhll(id)

                  enL = ((SL - qL)*UL(ie) - ptL*qL + pt*SM + &
                        Uhll(ibn)*(vBL - SM*Uhll(ibn) - at1*Uhll(ibt1) - at2*Uhll(ibt2))) &
                        / (SL - SM)
                  enR = ((SR - qR)*UR(ie) - ptR*qR + pt*SM + &
                        Uhll(ibn)*(vBR - SM*Uhll(ibn) - at1*Uhll(ibt1) - at2*Uhll(ibt2))) &
                        / (SR - SM)

                  ! weight factor, 0 or 1
                  f1 = 0.5e0 + sign(0.5e0, SM)
                  f2 = 1.0e0 - f1

                  Flux(id) =   f1*(SL*(dsL - UL(id)) + FL(id)) &
                           +   f2*(SR*(dsR - UR(id)) + FR(id))
                  Flux(imn) =  f1*(SL*(dsL*SM - UL(imn)) + FL(imn)) &
                            +  f2*(SR*(dsR*SM - UR(imn)) + FR(imn))
                  Flux(imt1) = f1*(SL*(dsL*at1 - UL(imt1)) + FL(imt1)) &
                             + f2*(SR*(dsR*at1 - UR(imt1)) + FR(imt1))
                  Flux(imt2) = f1*(SL*(dsL*at2 - UL(imt2)) + FL(imt2)) &
                             + f2*(SR*(dsR*at2 - UR(imt2)) + FR(imt2))
                  Flux(ie) =   f1*(SL*(enL - UL(ie)) + FL(ie)) &
                           +   f2*(SR*(enR - UR(ie)) + FR(ie))
                  else
                  ! HLLD flux
                  ! L* state
                  f1 = 1.e0 / (dL*(SL - qL)*(SL-SM) - Uhll(ibn)**2)
                  UL1(ibt1) = VLibt1 * f1 * (dL*(SL-qL)**2 - Uhll(ibn)**2)
                  UL1(ibt2) = VLibt2 * f1 * (dL*(SL-qL)**2 - Uhll(ibn)**2)
                  vt1L = VLimt1 - f1 * Uhll(ibn)*VLibt1*(SM-qL)
                  vt2L = VLimt2 - f1 * Uhll(ibn)*VLibt2*(SM-qL)

                  UL1(id) = dsL
                  UL1(ie) = ((SL-qL)*UL(ie) - ptL*qL + pt*SM + &
                    Uhll(ibn)*(vBL - SM*Uhll(ibn) - vt1L*UL1(ibt1) - vt2L*UL1(ibt2))) / (SL -SM)
                  UL1(imn) = dsL * SM
                  UL1(imt1) = dsL * vt1L
                  UL1(imt2) = dsL * vt2L

                  ! R* state
                  f1 = 1.e00 / (dR*(SR - qR)*(SR-SM) - Uhll(ibn)**2)
                  UR1(ibt1) = VRibt1 * f1 * (dR*(SR-qR)**2 - Uhll(ibn)**2)
                  UR1(ibt2) = VRibt2 * f1 * (dR*(SR-qR)**2 - Uhll(ibn)**2)
                  vt1R = VRimt1 - f1 * Uhll(ibn)*VRibt1*(SM-qR)
                  vt2R = VRimt2 - f1 * Uhll(ibn)*VRibt2*(SM-qR)

                  UR1(id) = dsR
                  UR1(ie) = ((SR-qR)*UR(ie) - ptR*qR + pt*SM + &
                    Uhll(ibn)*(vBR - SM*Uhll(ibn) - vt1R*UR1(ibt1) - vt2R*UR1(ibt2))) / (SR -SM)
                  UR1(imn) = dsR * SM
                  UR1(imt1) = dsR * vt1R
                  UR1(imt2) = dsR * vt2R

                  ! rotational waves
                  SL1 = SM - aVL
                  SR1 = SM + aVR

                  ! F = F(L*)
                  if (SL1 >= 0.000) then
                    Flux(id:ie) = SL * (UL1(id:ie) - UL(id:ie)) + FL(id:ie)
                    Flux(ibt1) =  SL * (UL1(ibt1) - UL(ibt1)) + FL(ibt1)
                    Flux(ibt2) =  SL * (UL1(ibt2) - UL(ibt2)) + FL(ibt2)

                  ! F = F(R*)
                  elseif (SR1 <= 0.000) then
                    Flux(id:ie) = SR * (UR1(id:ie) - UR(id:ie)) + FR(id:ie)
                    Flux(ibt1) =  SR * (UR1(ibt1) - UR(ibt1)) + FR(ibt1)
                    Flux(ibt2) =  SR * (UR1(ibt2) - UR(ibt2)) + FR(ibt2)

                  ! central states
                  else
                    sdsL = sqrt(dsL)
                    sdsR = sqrt(dsR)

                    f1 = 1.0e0/(sdsL + sdsR)
                    f2 = sign(1.0e0, Uhll(ibn))

                    at1 = f1 * (sdsL*vt1L + sdsR*vt1R + (UR1(ibt1) - UL1(ibt1))*f2)
                    at2 = f1 * (sdsL*vt2L + sdsR*vt2R + (UR1(ibt2) - UL1(ibt2))*f2)

                    U2(ibt1) = f1 * (sdsL * UR1(ibt1) + sdsR * UL1(ibt1) + sdsL*sdsR*(vt1R - vt1L)*f2)
                    U2(ibt2) = f1 * (sdsL * UR1(ibt2) + sdsR * UL1(ibt2) + sdsL*sdsR*(vt2R - vt2L)*f2)

                    ! F = F(L**)
                    if (SM >= 0.000) then

                      U2(id) = dsL
                      U2(ie) = UL1(ie) - sdsL*(vt1L*UL1(ibt1) + vt2L*UL1(ibt2) &
                                       - at1*U2(ibt1) - at2*U2(ibt2))*f2

                      U2(imn) = dsL * SM
                      U2(imt1) = dsL * at1
                      U2(imt2) = dsL * at2

                      Flux(id:ie) = SL1*U2(id:ie) - (SL1 - SL)*UL1(id:ie) - SL*UL(id:ie) + FL(id:ie)
                      !Flux(1:5) = SL1*U2(1:5) - (SL1 - SL)*UL1(1:5) - SL*UL(1:5) + FL(1:5)
                      Flux(ibt1) = SL1*U2(ibt1) - (SL1 - SL)*UL1(ibt1) - SL*UL(ibt1) + FL(ibt1)
                      Flux(ibt2) = SL1*U2(ibt2) - (SL1 - SL)*UL1(ibt2) - SL*UL(ibt2) + FL(ibt2)

                    ! F = F(R**)
                    else

                      U2(id) = dsR
                      U2(ie) = UR1(ie) + sdsR*(vt1R*UR1(ibt1) + vt2R*UR1(ibt2) &
                                       - at1*U2(ibt1) - at2*U2(ibt2))*f2

                      U2(imn) = dsR * SM
                      U2(imt1) = dsR * at1
                      U2(imt2) = dsR * at2

                      Flux(id:ie) = SR1*U2(id:ie) - (SR1 - SR)*UR1(id:ie) - SR*UR(id:ie) + FR(id:ie)
                      !Flux(1:5) = SL1*U2(1:5) - (SL1 - SL)*UL1(1:5) - SL*UL(1:5) + FL(1:5)
                      Flux(ibt1) = SR1*U2(ibt1) - (SR1 - SR)*UR1(ibt1) - SR*UR(ibt1) + FR(ibt1)
                      Flux(ibt2) = SR1*U2(ibt2) - (SR1 - SR)*UR1(ibt2) - SR*UR(ibt2) + FR(ibt2)

                    endif
                  endif
                  endif
              end if

              Flux(ibn) = 0.5e0 * ((VLips + VRips) - ch * (VRibn - VLibn))
              Flux(ips) = 0.5e0 * (ch*ch*(VLibn+VRibn) - ch*(VRips - VLips))

              if (dir == IAXIS) then
                xF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              else if (dir == JAXIS) then
                yF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              else if (dir == KAXIS) then
                zF(1:NCONSVAR_NUMBER,i,j,k) = Flux(1:NCONSVAR_NUMBER)
              end if
            end do
          end do
        end do
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        end do ! dir

        return
    end subroutine

#endif

    function eos_gete(U) result(ef)
        ! EOS implementation for the riemann solver
        use sim_data, only: gamma
        implicit none
        real, dimension(NCONSVAR_NUMBER) :: U
        real :: ef
        real :: df, uf, vf, wf, pf, bxf, byf, bzf

        df = U(DENS_VAR)
        uf = U(VELX_VAR)
        vf = U(VELY_VAR)
        wf = U(VELZ_VAR)
        pf = U(PRES_VAR)
        bxf = 0.0; byf = 0.0; bzf = 0.0
#ifdef MHD
        bxf = U(BMFX_VAR)
        byf = U(BMFY_VAR)
        bzf = U(BMFZ_VAR)
#endif
        ef = pf/(gamma - 1.0) + 0.5 * (df * (uf**2 + vf**2 + wf**2) + &
                                         (bxf**2 + byf**2 + bzf**2))
    end function eos_gete

end module flux_module

