module customBC_module
#include "param.h"
  implicit none
contains

  subroutine applycustomBC(var, dir, face)
                        
    use sim_data, only: nghost, ihi, ilo, jhi, jlo, khi, klo, &
                        jGlo, jGhi, iGlo, iGhi, kGlo, kGhi, PI, &
                        myrank, mainVar
    use grid_func, only: get_coords
    implicit none
    real(8), pointer :: q(:, :, :)
    integer, intent(in) :: var, Dir, Face
    real(8) :: sig
    integer :: ii, j, i, k
    real(8), allocatable, dimension(:) :: x,y,z

    !allocate(x(iGlo:iGhi), y(jGlo:jGhi), z(kGlo:kGhi))
    ! 
    !q(iGlo:, jGlo:, kGlo:) => mainVar(var,:,:,:)

    !call get_coords(IAXIS,iGlo,iGhi,x)
    !call get_coords(JAXIS,jGlo,jGhi,y)
    !call get_coords(KAXIS,kGlo,kGhi,z)

    !if ((Dir == IAXIS) .and. (Face == LEFT)) then
    !    ! Apply boundary at XL
    !end if
    !
    !if ((Dir == KAXIS) .and. (Face == LEFT)) then
    !    ! Apply boundary at ZL
    !    do ii = 1, nghost(KAXIS)
    !       q(:,:,klo-ii) = q(:,:,klo)
    !    end do
    !end if
    ! 
    !nullify(q) 
    !deallocate(x,y,z)
  end subroutine applycustomBC

end module customBC_module
