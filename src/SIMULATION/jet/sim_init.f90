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
    error stop "jet initialization currently supports 2D only"
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
       print *, "A Jet problem with custom X-left boundary selected"
    end if
    do k = klo, khi
      do j = jlo, jhi
        do i = ilo, ihi
           ! initialize the fields
           solnVar(DENS_VAR,i,j,k) = 1.0
           solnVar(PRES_VAR,i,j,k) = 2.5
           solnVar(VELX_VAR,i,j,k) = 0.0
           solnVar(VELY_VAR,i,j,k) = 0.0
           call eos_gete(solnVar(:,i,j,k))
        end do
      end do
    end do
    deallocate(x,y)
    nullify(solnVar)
  end subroutine init_problem
end module sim_init
