module applyBC_module
#include "param.h"
  implicit none
contains

  subroutine applydefBC(var, Dir, Face, bc_int)
    use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                        iGlo, jGlo, kGlo, mainVar, nghost
    implicit none
    real(8), pointer :: q(:,:,:)
    integer, intent(in) :: var, Dir, Face, bc_int
    integer :: lo(MAXDIM), hi(MAXDIM), ii, src, dest
    real(8) :: sig

    if (bc_int == CUSTOM .or. bc_int == PERIODIC) return
    if (Dir < IAXIS .or. Dir > MAXDIM) error stop "Invalid boundary direction"
    if (Dir > NDIM) return
    if (Face /= LEFT .and. Face /= RIGHT) error stop "Invalid boundary face"
    if (bc_int /= FLOW .and. bc_int /= REFLECT) error stop "Invalid boundary type"

    q(iGlo:,jGlo:,kGlo:) => mainVar(var,:,:,:)
    ! interior bounds in each direction
    lo = [ilo,jlo,klo]
    hi = [ihi,jhi,khi]

    ! reverse the normal velocity at reflecting walls
    sig = 1.0d0
    if (bc_int == REFLECT .and. var == VELX_VAR+Dir-1) sig = -1.0d0

    do ii = 1, nghost(Dir)
      if (Face == LEFT) then
        ! lower face: copy the edge cell for flow, mirror for reflect
        dest = lo(Dir) - ii
        src = lo(Dir)
        if (bc_int == REFLECT) src = lo(Dir) + ii - 1
      else
        ! upper face: copy the edge cell for flow, mirror for reflect
        dest = hi(Dir) + ii
        src = hi(Dir)
        if (bc_int == REFLECT) src = hi(Dir) - ii + 1
      end if

      ! fill the boundary plane
      select case (Dir)
      case (IAXIS)
        q(dest,:,:) = sig * q(src,:,:)
      case (JAXIS)
        q(:,dest,:) = sig * q(:,src,:)
      case (KAXIS)
        q(:,:,dest) = sig * q(:,:,src)
      end select
    end do

    nullify(q)
  end subroutine applydefBC

  subroutine applyBC_all()
    use sim_data, only: at_xlboundary, at_xrboundary, &
                        at_ylboundary, at_yrboundary, &
                        at_zlboundary, at_zrboundary, &
                        xlbc_int, xrbc_int, ylbc_int, yrbc_int, &
                        zlbc_int, zrbc_int
    use customBC_module, only: applycustomBC
    implicit none
    logical :: at_boundary(2,MAXDIM)
    integer :: bc_int(2,MAXDIM), dir, side, face, n

    ! Periodic BC is already applied in the guardcell_fill routine
    ! applyBC_all should always be called after guardcell_fill routine

    ! lower and upper faces in each direction
    at_boundary(:,IAXIS) = [at_xlboundary,at_xrboundary]
    at_boundary(:,JAXIS) = [at_ylboundary,at_yrboundary]
    at_boundary(:,KAXIS) = [at_zlboundary,at_zrboundary]
    bc_int(:,IAXIS) = [xlbc_int,xrbc_int]
    bc_int(:,JAXIS) = [ylbc_int,yrbc_int]
    bc_int(:,KAXIS) = [zlbc_int,zrbc_int]

    ! apply boundaries only in active directions
    do dir = 1, NDIM
      do side = 1, 2
        if (.not. at_boundary(side,dir)) cycle
        ! side 1 is LEFT, side 2 is RIGHT
        face = merge(LEFT,RIGHT,side == 1)
        do n = NVAR_BEGIN, NVAR_END
          if (bc_int(side,dir) == CUSTOM) then
            call applycustomBC(n,dir,face)
          else
            call applydefBC(n,dir,face,bc_int(side,dir))
          end if
        end do
      end do
    end do
  end subroutine applyBC_all
end module applyBC_module
