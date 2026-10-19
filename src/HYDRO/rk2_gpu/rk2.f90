module rk2_module
#include "param.h"
  implicit none
contains

    subroutine RK2_SSP(dt)
       use flux_module, only: riemann_solve
       use recon_module, only: recon_getcellfaces
       use eos_module, only: eos_getp
       use sim_data, only: grav, usegrav, ilo, ihi, jlo, jhi, klo, khi, &
                           dx, dy, dz, mainVar, iGlo, jGlo, kGlo, &
                           iGhi, jGhi, kGhi, smalld, smallp, smalle, gamma

       use applyBC_module, only: applyBC_all
       use guard_func, only: guardcell_fill
       implicit none
       real(8), pointer :: solnVar(:,:,:,:)
       real(8), intent(in) :: dt

       real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi) :: U, Up1
       real(8), dimension(1:NDIM,1:NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi) :: recon_plus, recon_minus

       real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi,kGlo:kGhi) :: xF, yF, zF
       real(8), dimension(1:NCONSVAR_NUMBER) :: sgrav, sterm

       integer :: i, j, k, stage, n

       solnVar(1:,iGlo:,jGlo:,kGlo:) => mainVar(1:,:,:,:)

       sterm(:) = 0.0
       if (usegrav) then
         sterm(3) = grav
         sterm(5) = grav
       end if

       !$OMP TARGET DATA MAP(ALWAYS, TOFROM: solnVar) &
       !$OMP MAP(TO: dt, gamma, smalld, smallp, smalle, sterm, dx, dy, dz) &
       !$OMP MAP(ALLOC: recon_plus, recon_minus, xF, yF, zF, U, Up1)

       ! Initialize the conserved arrays
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)
       do k = kGlo, kGhi
         do j = jGlo, jGhi
           do i = iGlo, iGhi
             Up1(DENS_CONS,i,j,k) = max(solnVar(DENS_VAR,i,j,k), smalld)
             Up1(MOMX_CONS:MOMZ_CONS,i,j,k) = solnVar(DENS_VAR,i,j,k) * solnVar(VELX_VAR:VELZ_VAR,i,j,k)
             Up1(ENER_CONS,i,j,k) = max(solnVar(ENER_VAR,i,j,k), smalle)

#ifdef MHD
             Up1(BMFX_CONS:BMFZ_CONS,i,j,k) = solnVar(BMFX_VAR:BMFZ_VAR,i,j,k)
             Up1(BPSI_CONS,i,j,k) = solnVar(BPSI_VAR,i,j,k)
#endif
             U(DENS_CONS,i,j,k) = max(solnVar(DENS_VAR,i,j,k), smalld)
             U(MOMX_CONS:MOMZ_CONS,i,j,k) = solnVar(DENS_VAR,i,j,k) * solnVar(VELX_VAR:VELZ_VAR,i,j,k)
             U(ENER_CONS,i,j,k) = max(solnVar(ENER_VAR,i,j,k), smalle)

#ifdef MHD
             U(BMFX_CONS:BMFZ_CONS,i,j,k) = solnVar(BMFX_VAR:BMFZ_VAR,i,j,k)
             U(BPSI_CONS,i,j,k) = solnVar(BPSI_VAR,i,j,k)
#endif
           end do
         end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

       ! rk2 has 2 steps
       do stage = 1, 2
         ! use the updated solution for the next step

         call recon_getcellfaces(solnVar, &
                                 recon_plus, recon_minus)
         call riemann_solve(recon_plus, recon_minus, xF, yF, zF)


          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP PRIVATE(sgrav)
          do k = klo, khi
            do j = jlo, jhi
              do i = ilo, ihi
#ifdef MHD
                sgrav = sterm * (/1.0e0, 0.0, Up1(DENS_CONS,i,j,k), 0.0, Up1(MOMY_CONS,i,j,k), &
                                  0.0,0.0,0.0,0.0/)
#else
                sgrav = sterm * (/1.0e0, 0.0, Up1(DENS_CONS,i,j,k), 0.0, Up1(MOMY_CONS,i,j,k)/)
#endif
                Up1(:,i,j,k) = Up1(:,i,j,k) &
                             + (dt/dx) * (xF(:,i,j,k) - xF(:,i+1,j,k)) &
#if NDIM >= 2
                             + (dt/dy) * (yF(:,i,j,k) - yF(:,i,j+1,k)) &
#endif
#if NDIM == 3
                             + (dt/dz) * (zF(:,i,j,k) - zF(:,i,j,k+1)) &
#endif
                             + dt * sgrav(:)

                ! get primitive quantities
                solnVar(DENS_VAR,i,j,k) = max(Up1(DENS_CONS,i,j,k), smalld)
                solnVar(VELX_VAR:VELZ_VAR,i,j,k) = Up1(MOMX_CONS:MOMZ_CONS,i,j,k)/Up1(DENS_CONS,i,j,k)
                solnVar(ENER_VAR,i,j,k) = max(Up1(ENER_CONS,i,j,k), smalle)
#ifdef MHD
                solnVar(BMFX_VAR:BMFZ_VAR,i,j,k) = Up1(BMFX_CONS:BMFZ_CONS,i,j,k)
                solnVar(BPSI_VAR,i,j,k) = Up1(BPSI_CONS,i,j,k)
#endif
                ! Apply EOS: in 2D
                ! the same expression applies in 1D and 3D
                solnVar(PRES_VAR,i,j,k) = (gamma - 1.0) * (solnVar(ENER_VAR,i,j,k) - 0.50 * solnVar(DENS_VAR,i,j,k) * &
                                        (solnVar(VELX_VAR,i,j,k)**2 + solnVar(VELY_VAR,i,j,k)**2 + &
                                         solnVar(VELZ_VAR,i,j,k)**2))
#ifdef MHD
                solnVar(PRES_VAR,i,j,k) = solnVar(PRES_VAR,i,j,k) - (gamma - 1.0) * 0.50 * &
                                        (solnVar(BMFX_VAR,i,j,k)**2 + solnVar(BMFY_VAR,i,j,k)**2 + &
                                         solnVar(BMFZ_VAR,i,j,k)**2)
#endif
                solnVar(PRES_VAR,i,j,k) = max(solnVar(PRES_VAR,i,j,k), smallp)
              end do
            end do
          end do
          !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
          !$OMP TARGET UPDATE FROM(solnVar)
          call guardcell_fill()
          call applyBC_all()
          !$OMP TARGET UPDATE TO(solnVar)

       end do ! rk step loop

       ! get the updated solution and store it as Up1

       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)
       do k = kGlo, kGhi
         do j = jGlo, jGhi
           do i = iGlo, iGhi
             Up1(DENS_CONS,i,j,k) = max(solnVar(DENS_VAR,i,j,k), smalld)
             Up1(MOMX_CONS:MOMZ_CONS,i,j,k) = solnVar(DENS_VAR,i,j,k) * solnVar(VELX_VAR:VELZ_VAR,i,j,k)
             Up1(ENER_CONS,i,j,k) = max(solnVar(ENER_VAR,i,j,k), smalle)

#ifdef MHD
             Up1(BMFX_CONS:BMFZ_CONS,i,j,k) = solnVar(BMFX_VAR:BMFZ_VAR,i,j,k)
             Up1(BPSI_CONS,i,j,k) = solnVar(BPSI_VAR,i,j,k)
#endif
           end do
         end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)
       do k = kGlo, kGhi
         do j = jGlo, jGhi
           do i = iGlo, iGhi
             Up1(:,i,j,k) = 0.5e0 * (U(:,i,j,k) + Up1(:,i,j,k))

             ! get primitive quantities
             solnVar(DENS_VAR,i,j,k) = max(Up1(DENS_CONS,i,j,k), smalld)
             solnVar(VELX_VAR:VELZ_VAR,i,j,k) = Up1(MOMX_CONS:MOMZ_CONS,i,j,k)/Up1(DENS_CONS,i,j,k)
             solnVar(ENER_VAR,i,j,k) = max(Up1(ENER_CONS,i,j,k), smalle)
#ifdef MHD
             solnVar(BMFX_VAR:BMFZ_VAR,i,j,k) = Up1(BMFX_CONS:BMFZ_CONS,i,j,k)
             solnVar(BPSI_VAR,i,j,k) = Up1(BPSI_CONS,i,j,k)
#endif
             ! Apply EOS: in 2D
             ! the same expression applies in 1D and 3D
             solnVar(PRES_VAR,i,j,k) = (gamma - 1.0) * (solnVar(ENER_VAR,i,j,k) - 0.50 * solnVar(DENS_VAR,i,j,k) * &
                                     (solnVar(VELX_VAR,i,j,k)**2 + solnVar(VELY_VAR,i,j,k)**2 + &
                                      solnVar(VELZ_VAR,i,j,k)**2))
#ifdef MHD
             solnVar(PRES_VAR,i,j,k) = solnVar(PRES_VAR,i,j,k) - (gamma - 1.0) * 0.50 * &
                                     (solnVar(BMFX_VAR,i,j,k)**2 + solnVar(BMFY_VAR,i,j,k)**2 + &
                                      solnVar(BMFZ_VAR,i,j,k)**2)
#endif
             solnVar(PRES_VAR,i,j,k) = max(solnVar(PRES_VAR,i,j,k), smallp)

           end do
         end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
       !$OMP END TARGET DATA

       ! refresh guard cells after the final RK combination
       call guardcell_fill()
       call applyBC_all()

       nullify(solnVar)
    end subroutine RK2_SSP
end module rk2_module
