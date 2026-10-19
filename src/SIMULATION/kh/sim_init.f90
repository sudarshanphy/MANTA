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
    integer :: i, j, k
    real, allocatable :: x(:), y(:)

#if NDIM != 2
    error stop "kh initialization currently supports 2D only"
#endif

    allocate(x(ilo:ihi))
    allocate(y(jlo:jhi))

    x = 0.0
    y = 0.0

    call get_coords(IAXIS,ilo,ihi,x)
    call get_coords(JAXIS,jlo,jhi,y)

    solnVar(1:,ilo:,jlo:,klo:) => mainVar(1:,ilo:ihi,jlo:jhi,klo:khi)
    solnVar(:,:,:,:) = 0.0

    if (myrank == MASTER_PROC) then
       print *, "Kelvin-Helmholtz (KH) problem selected"
    end if
    do k = klo, khi
      do j = jlo, jhi
        do i = ilo, ihi
           ! initialize the fields
           if (abs(y(j) - 0.5) > 0.25) then
             solnVar(DENS_VAR, i,j,k) = 1.0e0
             solnVar(VELX_VAR, i,j,k) = -0.5e0
             solnVar(VELY_VAR, i,j,k) = 0.0e0
             solnVar(PRES_VAR, i,j,k) = 2.5e0
           else
             solnVar(DENS_VAR, i,j,k) = 2.0e0
             solnVar(VELX_VAR, i,j,k) = 0.5e0
             solnVar(VELY_VAR, i,j,k) = 0.0e0
             solnVar(PRES_VAR, i,j,k) = 2.5e0
           end if
           ! give velocity perturbations
           if (abs(y(j) - 0.25) < 0.1) then
              solnVar(VELY_VAR, i,j,k) = 0.05 * sin(2.0 * 2.0 * PI * x(i))
           end if
           if (abs(y(j) - 0.75) < 0.1) then
              solnVar(VELY_VAR, i,j,k) = 0.05 * sin(2.0 * 2.0 * PI * x(i))
           end if
           call eos_gete(solnVar(:,i,j,k))
        end do
      end do
    end do
    deallocate(x,y)
    nullify(solnVar)
  end subroutine init_problem
end module sim_init
