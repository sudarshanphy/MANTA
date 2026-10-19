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
    real, allocatable :: x(:), y(:), z(:)

    allocate(x(ilo:ihi))
    allocate(y(jlo:jhi))
    allocate(z(klo:khi))

    x = 0.0
    y = 0.0
    z = 0.0

    call get_coords(IAXIS,ilo,ihi,x)
    call get_coords(JAXIS,jlo,jhi,y)
    call get_coords(KAXIS,klo,khi,z)

    solnVar(1:,ilo:,jlo:,klo:) => mainVar(1:,ilo:ihi,jlo:jhi,klo:khi)
    solnVar(:,:,:,:) = 0.0

    if (myrank == MASTER_PROC) then
       print *, "Sedov problem selected"
    end if
    do k = klo, khi
      do j = jlo, jhi
        do i = ilo, ihi
           ! initialize the fields
           solnVar(DENS_VAR, i,j,k) = 1.0e0
           solnVar(VELX_VAR, i,j,k) = 0.0e0
           solnVar(VELY_VAR, i,j,k) = 0.0e0
           solnVar(PRES_VAR, i,j,k) = 1.0e-1

           ! add a region with very high pressure
#if NDIM == 1
           r = abs(x(i))
#elif NDIM == 2
           r = sqrt(x(i)**2 + y(j)**2)
#else
           r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
#endif
           if (r <= 0.1) then
             solnVar(PRES_VAR, i,j,k) = 1.0e1
           end if
#ifdef MHD
           solnVar(BMFX_VAR, i,j,k) = 1.0/sqrt(2.0)
           solnVar(BMFY_VAR, i,j,k) = 1.0/sqrt(2.0)
           solnVar(BMFZ_VAR, i,j,k) = 0.0
#endif
           ! compute energy after all fields are initialized
           call eos_gete(solnVar(:,i,j,k))
        end do
      end do
    end do
    deallocate(x,y,z)
    nullify(solnVar)
  end subroutine init_problem
end module sim_init
