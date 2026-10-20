module flux_module
#include "param.h"
    implicit none

    ! HLLC solver is only for pure hydro simulation
    ! HLLE can be used for both pure Hydro and MHD simulation
    ! HLLD solver is only for pure MHD simulation
contains

#ifndef MHD
    subroutine hllc(Uleft, Uright, Dir, Flux)
        use sim_data, only: gamma
        use misc_module, only: to_upper 
        implicit none

        real(8), dimension(1:NCONSVAR_NUMBER), intent(in) :: Uleft, Uright
        character(len=1), intent(in) :: Dir
        real(8), dimension(1:NCONSVAR_NUMBER), intent(out) :: Flux
        real(8), dimension(1:NCONSVAR_NUMBER) :: FL, FR, UL, UR, FstarL, FstarR, UstarL, UstarR
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

        dL = Uleft(1)
        ufL = Uleft(2)
        vfL = Uleft(3)
        wfL = Uleft(4)
        pL = Uleft(5)

        dR = Uright(1)
        ufR = Uright(2)
        vfR = Uright(3)
        wfR = Uright(4)
        pR = Uright(5)

        select case (to_upper(Dir))
        case ('X')
            n = (/1, 0, 0/)
        case ('Y')
            n = (/0, 1, 0/)
        case ('Z')
            n = (/0, 0, 1/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

        sdR = sqrt(dR)
        sdL = sqrt(dL)

        udL = ufL * dL
        udR = ufR * dR

        vdL = vfL * dL
        vdR = vfR * dR

        wdL = wfL * dL
        wdR = wfR * dR

        eL = eos_gete((/dL, ufL, vfL, wfL, pL/))
        eR = eos_gete((/dR, ufR, vfR, wfR, pR/))

        HL = (eL + pL) / dL
        HR = (eR + pR) / dR

        cL = sqrt(gamma * pL / dL)
        cR = sqrt(gamma * pR / dR)

        ubar = (ufL * sdL + ufR * sdR) / (sdL + sdR)
        vbar = (vfL * sdL + vfR * sdR) / (sdL + sdR)
        wbar = (wfL * sdL + wfR * sdR) / (sdL + sdR)

        Hbar = (HL * sdL + HR * sdR) / (sdL + sdR)

        cbar = sqrt((gamma - 1.00) * (Hbar - 0.50 * (ubar**2 + vbar**2 + wbar**2)))

        qL = sum((/ufL, vfL, wfL/) * n)
        qR = sum((/ufR, vfR, wfR/) * n)

        qbar = sum((/ubar, vbar, wbar/) * n)

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
        return
    end subroutine hllc
#endif

    subroutine hlle(Uleft, Uright, Dir, Flux)
        use sim_data, only: gamma
#ifdef MHD
        use sim_data, only: ch
#endif
        use misc_module, only: to_upper 
        implicit none

        character(len=1), intent(in) :: Dir
        real(8), dimension(1:NCONSVAR_NUMBER), intent(in) :: Uleft, Uright
        real(8), dimension(1:NCONSVAR_NUMBER), intent(out) :: Flux
        real(8), dimension(1:NCONSVAR_NUMBER) :: FL, FR, UL, UR, Fstar
        real(8) :: dL,  pL 
        real(8) :: dR,  pR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, HL, HR, cL, cR
#ifdef MHD
        real(8) :: bxL, byL, bzL, bpL
        real(8) :: bxR, byR, bzR, bpR
        real(8) :: B2L, B2R, vBL, vBR, qBL, qBR, cfL, cfR
        real(8) :: bnstar, psistar
#endif
        real(8) :: ubar, vbar, wbar, Hbar, cbar
        real(8) :: qL, qR, qbar, SL, SR
        real(8), dimension(3) :: n
        ! index for different dimension
        integer :: id, imn, imt1, imt2, ibn, ibt1, ibt2, ie, ips

        dL = Uleft(1)
        ufL = Uleft(2)
        vfL = Uleft(3)
        wfL = Uleft(4)
        pL = Uleft(5)

        dR = Uright(1)
        ufR = Uright(2)
        vfR = Uright(3)
        wfR = Uright(4)
        pR = Uright(5)

#ifdef MHD
        bxL = Uleft(6)
        byL = Uleft(7)
        bzL = Uleft(8)
        bpL = Uleft(9)
        
        bxR = Uright(6)
        byR = Uright(7)
        bzR = Uright(8)
        bpR = Uright(9)
#endif
        id = 1; ie = 5; ips = 9
        select case (to_upper(Dir))
        case ('X')
            imn = 2; imt1 = 3; imt2 = 4
            ibn = 6; ibt1 = 7; ibt2 = 8
            n = (/1, 0, 0/)
        case ('Y')
            imn = 3; imt1 = 2; imt2 = 4
            ibn = 7; ibt1 = 6; ibt2 = 8
            n = (/0, 1, 0/)
        case ('Z')
            imn = 4; imt1 = 2; imt2 = 3
            ibn = 8; ibt1 = 6; ibt2 = 7
            n = (/0, 0, 1/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

#ifdef MHD
        ! solve the cleaning waves before the MHD Riemann problem
        ! GLM interface construction: Mignone & Tzeferacos, JCP 229 (2010), 2117-2138.
        ! "A Second-Order Unsplit Godunov Scheme for Cell-Centered MHD: the CTU-GLM scheme"
        ! Sec. 3.2, Eqs. (23)-(25): https://arxiv.org/html/0911.3410#S3.SS2
        ! Only this interface construction is used here, with our HLLE/HLLD solver.
        ! Bn is the face-normal field; psi is the cleaning variable (bpL, bpR).
        ! d_t Bn + d_n psi = 0,   d_t psi + ch**2 d_n Bn = 0.
        ! Bn*  = (BnL+BnR)/2 - (psiR-psiL)/(2*ch)
        ! psi* = (psiL+psiR)/2 - ch*(BnR-BnL)/2
        ! Previously only the Bn and psi fluxes were corrected after the MHD solve.
        ! Now both MHD states use Bn* and psi* before energies and wave speeds are built.
        ! Density, velocity, pressure and transverse magnetic fields are unchanged.
        bnstar = 0.5*(Uleft(ibn)+Uright(ibn)) - 0.5*(bpR-bpL)/ch
        psistar = 0.5*(bpL+bpR) - 0.5*ch*(Uright(ibn)-Uleft(ibn))
        select case (ibn)
        case (6)
          bxL = bnstar; bxR = bnstar
        case (7)
          byL = bnstar; byR = bnstar
        case (8)
          bzL = bnstar; bzR = bnstar
        end select
        bpL = psistar; bpR = psistar
#endif

        sdR = sqrt(dR)
        sdL = sqrt(dL)

        cL = sqrt(gamma * pL / dL)
        cR = sqrt(gamma * pR / dR)
        
        ! speed in the normal direction
        qL = Uleft(imn)
        qR = Uright(imn)

#ifdef MHD
        ! rebuild total energy with Bn* while keeping pressure fixed:
        ! E_s* = E_s + ((Bn*)**2-Bn_s**2)/2, s=L,R (before floors).
        eL =  eos_gete((/dL, ufL, vfL, wfL, pL, bxL, byL, bzL, bpL/))
        eR =  eos_gete((/dR, ufR, vfR, wfR, pR, bxR, byR, bzR, bpR/))
#else
        eL =  eos_gete(Uleft(1:5))
        eR =  eos_gete(Uright(1:5))
#endif

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
        qBL = bnstar
        qBR = bnstar
        
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
        ! cleaning fluxes: F_Bn = psi*, F_psi = ch**2*Bn*.
        ! These equal the previous cleaning fluxes; the other MHD fluxes now use Bn*.
        Flux(ips) = ch*ch*bnstar
        Flux(ibn) = psistar
#endif
        return
    end subroutine hlle

#ifdef MHD
    subroutine hlld(VL, VR, Dir, Flux)
        use sim_data, only: gamma, ch
        use misc_module, only: to_upper 
        implicit none

        real(8), dimension(1:NCONSVAR_NUMBER), intent(in) :: VL, VR
        character(len=1), intent(in) :: Dir
        real(8), dimension(1:NCONSVAR_NUMBER), intent(out) :: Flux
        real(8), dimension(1:NCONSVAR_NUMBER) :: FL, FR, UL, UR, UL1, U2, UR1, Uhll
        real(8) :: dL, udL, vdL, wdL, pL, bxL, byL, bzL, bpL
        real(8) :: dR, udR, vdR, wdR, pR, bxR, byR, bzR, bpR
        real(8) :: ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, cfL, cfR, ptL, ptR, pt 
        real(8) :: qL, qR, SL, SR, SM, qBL, qBR, B2L, B2R, vBL, vBR
        real(8) :: bnstar, psistar
        real(8) :: f1, f2, SL1, SR1, aVL, aVR, at1, at2, vt1L, vt2L
        real(8) :: vt1R, vt2R, dsL, dsR, sdsL, sdsR, enL, enR
        real(8), dimension(3) :: n
        real(8), parameter :: hllg_factor = 1.001d0
        ! index for different dimension
        integer :: id, imn, imt1, imt2, ibn, ibt1, ibt2, ie, ips

        dL =  VL(1)
        ufL = VL(2)
        vfL = VL(3)
        wfL = VL(4)
        pL =  VL(5)
        bxL = VL(6)
        byL = VL(7)
        bzL = VL(8)
        bpL = VL(9)

        dR =  VR(1)
        ufR = VR(2)
        vfR = VR(3)
        wfR = VR(4)
        pR =  VR(5)
        bxR = VR(6)
        byR = VR(7)
        bzR = VR(8)
        bpR = VR(9)

        id = 1; ie = 5; ips = 9
        select case (to_upper(Dir))
        case ('X')
            imn = 2; imt1 = 3; imt2 = 4
            ibn = 6; ibt1 = 7; ibt2 = 8
            n = (/1.0, 0.0, 0.0/)
        case ('Y')
            imn = 3; imt1 = 2; imt2 = 4
            ibn = 7; ibt1 = 6; ibt2 = 8
            n = (/0.0, 1.0, 0.0/)
        case ('Z')
            imn = 4; imt1 = 2; imt2 = 3
            ibn = 8; ibt1 = 6; ibt2 = 7
            n = (/0.0, 0.0, 1.0/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

#ifdef MHD
        ! solve the cleaning waves before the MHD Riemann problem
        ! GLM interface construction: Mignone & Tzeferacos, JCP 229 (2010), 2117-2138.
        ! "A Second-Order Unsplit Godunov Scheme for Cell-Centered MHD: the CTU-GLM scheme"
        ! Sec. 3.2, Eqs. (23)-(25): https://arxiv.org/html/0911.3410#S3.SS2
        ! Only this interface construction is used here, with our HLLE/HLLD solver.
        ! Bn is the face-normal field; psi is the cleaning variable (bpL, bpR).
        ! d_t Bn + d_n psi = 0,   d_t psi + ch**2 d_n Bn = 0.
        ! Bn*  = (BnL+BnR)/2 - (psiR-psiL)/(2*ch)
        ! psi* = (psiL+psiR)/2 - ch*(BnR-BnL)/2
        ! Previously only the Bn and psi fluxes were corrected after the MHD solve.
        ! Now both MHD states use Bn* and psi* before energies and wave speeds are built.
        ! Density, velocity, pressure and transverse magnetic fields are unchanged.
        bnstar = 0.5*(VL(ibn)+VR(ibn)) - 0.5*(bpR-bpL)/ch
        psistar = 0.5*(bpL+bpR) - 0.5*ch*(VR(ibn)-VL(ibn))
        select case (ibn)
        case (6)
          bxL = bnstar; bxR = bnstar
        case (7)
          byL = bnstar; byR = bnstar
        case (8)
          bzL = bnstar; bzR = bnstar
        end select
        bpL = psistar; bpR = psistar
#endif

        udL = ufL * dL
        udR = ufR * dR

        vdL = vfL * dL
        vdR = vfR * dR

        wdL = wfL * dL
        wdR = wfR * dR

        ! rebuild total energy with Bn* while keeping pressure fixed:
        ! E_s* = E_s + ((Bn*)**2-Bn_s**2)/2, s=L,R (before floors).
        eL =  eos_gete((/dL, ufL, vfL, wfL, pL, bxL, byL, bzL, bpL/))
        eR =  eos_gete((/dR, ufR, vfR, wfR, pR, bxR, byR, bzR, bpR/))

        qL = VL(imn) 
        qR = VR(imn)

        ! magnetic field value in the normal direction
        qBL = bnstar
        qBR = bnstar
        
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
         
        !print *, "Inside HLLD" 

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
           pt = ptL + VL(id) * (SL - VL(imn)) * (SM - VL(imn))
           ! density in star region: M05 eqn 43
           dsL = VL(id) * (SL - VL(imn)) / (SL - SM) !d(L*)
           dsR = VR(id) * (SR - VR(imn)) / (SR - SM) !d(R*)
           
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

         if ((SL >= (SM - hllg_factor * aVL)) .or. &
               ((SM + hllg_factor * aVR) >= SR)) then
             f1 = 1.0e0 / (SR - SL)
             f2 = SR * SL
             Flux(ibt1) = f1 * (SR * FL(ibt1) - SL * FR(ibt1) + f2 * (VR(ibt1) - VL(ibt1)))  
             Flux(ibt2) = f1 * (SR * FL(ibt2) - SL * FR(ibt2) + f2 * (VR(ibt2) - VL(ibt2))) 

             at1 = Uhll(imt1)/Uhll(id)
             at2 = Uhll(imt2)/Uhll(id)
             
             enL = ((SL - VL(imn))*UL(ie) - ptL*VL(imn) + pt*SM + &
                   Uhll(ibn)*(vBL - SM*Uhll(ibn) - at1*Uhll(ibt1) - at2*Uhll(ibt2))) &
                   / (SL - SM)
             enR = ((SR - VR(imn))*UR(ie) - ptR*VR(imn) + pt*SM + &
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
           f1 = 1.e0 / (VL(id)*(SL - VL(imn))*(SL-SM) - Uhll(ibn)**2)
           UL1(ibt1) = VL(ibt1) * f1 * (VL(id)*(SL-VL(imn))**2 - Uhll(ibn)**2)
           UL1(ibt2) = VL(ibt2) * f1 * (VL(id)*(SL-VL(imn))**2 - Uhll(ibn)**2)
           vt1L = VL(imt1) - f1 * Uhll(ibn)*VL(ibt1)*(SM-VL(imn))
           vt2L = VL(imt2) - f1 * Uhll(ibn)*VL(ibt2)*(SM-VL(imn))
             
           UL1(id) = dsL
           UL1(ie) = ((SL-VL(imn))*UL(ie) - ptL*VL(imn) + pt*SM + &
             Uhll(ibn)*(vBL - SM*Uhll(ibn) - vt1L*UL1(ibt1) - vt2L*UL1(ibt2))) / (SL -SM)
           UL1(imn) = dsL * SM
           UL1(imt1) = dsL * vt1L
           UL1(imt2) = dsL * vt2L

           ! R* state
           f1 = 1.e00 / (VR(id)*(SR - VR(imn))*(SR-SM) - Uhll(ibn)**2)
           UR1(ibt1) = VR(ibt1) * f1 * (VR(id)*(SR-VR(imn))**2 - Uhll(ibn)**2)
           UR1(ibt2) = VR(ibt2) * f1 * (VR(id)*(SR-VR(imn))**2 - Uhll(ibn)**2)
           vt1R = VR(imt1) - f1 * Uhll(ibn)*VR(ibt1)*(SM-VR(imn))
           vt2R = VR(imt2) - f1 * Uhll(ibn)*VR(ibt2)*(SM-VR(imn))
             
           UR1(id) = dsR
           UR1(ie) = ((SR-VR(imn))*UR(ie) - ptR*VR(imn) + pt*SM + &
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
               U2(ie) = UR1(ie) - sdsR*(vt1R*UR1(ibt1) + vt2R*UR1(ibt2) &
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

        ! cleaning fluxes: F_Bn = psi*, F_psi = ch**2*Bn*.
        ! These equal the previous cleaning fluxes; the other MHD fluxes now use Bn*.
        Flux(ibn) = psistar
        Flux(ips) = ch*ch*bnstar

        return
    end subroutine hlld
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

