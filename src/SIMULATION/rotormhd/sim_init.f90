module sim_init
#include "param.h"
  implicit none
contains
  subroutine init_problem()

    use sim_data, only: ihi, ilo, jhi, &
                        jlo, klo, khi, PI, grav, usegrav, gamma, &
                        mainVar, myrank
    use grid_func, only: get_coords
    use eos_module, only: eos_gete, eos_getp

    implicit none
    real, pointer :: solnVar(:,:,:,:)
    integer :: i, j, k, iInt, jInt
    real :: r
    real, allocatable :: x(:), y(:)

#if NDIM != 2
    error stop "rotormhd initialization currently supports 2D only"
#endif

    allocate(x(ilo:ihi))
    allocate(y(jlo:jhi))

    x = 0.0
    y = 0.0

    call get_coords(IAXIS,ilo,ihi,x)
    call get_coords(JAXIS,jlo,jhi,y)

    solnVar(1:,ilo:,jlo:,klo:) => mainVar(1:,ilo:ihi,jlo:jhi,klo:khi)
    solnVar(:,:,:,:) = 0.0

    ! setup take from "Limiting and divergence cleaning for continuous finite
    ! element discretizations of the MHD equations" by Dmitri Kuzmin
    ! page no 18
    if (myrank == MASTER_PROC) then
       print *, "2D ROTOR MHD Problem"
    end if
    do k = klo, khi
      do j = jlo, jhi
        do i = ilo, ihi
           r = sqrt((x(i) - 0.5)**2 + (y(j) - 0.5)**2)
           if (r < 0.1) then
             solnVar(DENS_VAR, i,j,k) = 10.0
             solnVar(VELX_VAR, i,j,k) = 10.0 * (0.5 - y(j))
             solnVar(VELY_VAR, i,j,k) = 10.0 * (x(i) - 0.5)
           else if (r > 0.115) then
             solnVar(DENS_VAR, i,j,k) = 1.0
             solnVar(VELX_VAR, i,j,k) = 0.0
             solnVar(VELY_VAR, i,j,k) = 0.0
           else
             solnVar(DENS_VAR, i,j,k) = 1.0 + 9.0 * (0.115 - r)/(r - 0.1)
             solnVar(VELX_VAR, i,j,k) =  100.0 * (0.115 - r)/(r - 0.1) * (0.5 - y(j)) &
                                       / solnVar(DENS_VAR, i,j,k)
             solnVar(VELY_VAR, i,j,k) =  100.0 * (0.115 - r)/(r - 0.1) * (x(i) - 0.5) &
                                       / solnVar(DENS_VAR, i,j,k)
             ! need this to make sure density doen't blow up at initialization
             solnVar(DENS_VAR, i,j,k) = min(solnVar(DENS_VAR, i,j,k), 10.0)
           end if
           solnVar(PRES_VAR, i,j,k) = 0.5
#ifdef MHD
           solnVar(BMFX_VAR, i,j,k) = 2.5 / sqrt(4.0 * PI)
           solnVar(BMFY_VAR, i,j,k) = 0.0
#endif
           call eos_gete(solnVar(:,i,j,k))
        end do
      end do
    end do
    deallocate(x,y)
    nullify(solnVar)
  end subroutine init_problem
end module sim_init
