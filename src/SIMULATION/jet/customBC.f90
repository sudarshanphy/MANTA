module customBC_module
#include "param.h"
  implicit none
contains

  subroutine applycustomBC(var, dir, face)
                        
    use sim_data, only: Gpts, ihi, ilo, jhi, jlo, &
                        jGlo, jGhi, iGlo, iGhi, kGlo, kGhi, PI, &
                        myrank, mainVar
    use grid_func, only: get_coords
    use eos_module, only: eos_gete
    implicit none
    real, pointer :: q(:, :, :)
    integer, intent(in) :: var, Dir, Face
    real :: sig
    integer :: ii, j, i, k
    real, allocatable, dimension(:) :: x,y

    allocate(x(iGlo:iGhi), y(jGlo:jGhi))
     
    q(iGlo:, jGlo:, kGlo:) => mainVar(var,:,:,:)

    call get_coords(IAXIS,iGlo,iGhi,x)
    call get_coords(JAXIS,jGlo,jGhi,y)

    if ((Dir == IAXIS) .and. (Face == LEFT)) then
        ! Apply boundary at XL
       sig = 1.000
       if (var == VELX_VAR) sig = -1.000
        ! reflecting boundary
       do ii = 1, Gpts
         ! lower face
         q(ilo-ii, :, :) = sig * q(ilo+ii-1, :, :)
       end do
       ! inside a certain region
       do j = jGlo, jGhi
         do i = iGlo, iGhi
           if (abs(y(j)) <= 0.01) then
             if (x(i) <= 0.25) then
                if (var==VELX_VAR) q(i,j,:) = 0.5
                if (var==DENS_VAR) q(i,j,:) = 2.0
                ! update energy after the inlet fields have been set
                if (var==ENER_VAR) then
                  do k = kGlo, kGhi
                    call eos_gete(mainVar(:,i,j,k))
                  end do
                end if
             end if
           end if
         end do
       end do
    end if
     
    nullify(q) 
    deallocate(x,y)
  end subroutine applycustomBC

end module customBC_module
