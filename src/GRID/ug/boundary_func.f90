module boundary_func
#include "param.h"

  implicit none
  contains

  subroutine boundary_init(ilo,ihi,jlo,jhi,klo,khi,nx,ny,nz)
    use sim_data, only: at_xlboundary, at_xrboundary, &
                        at_ylboundary, at_yrboundary, &
                        at_zlboundary, at_zrboundary, &
                        xrank, yrank, zrank, xlbctype, xrbctype, &
                        ylbctype, yrbctype, xlbc_int, xrbc_int, &
                        ylbc_int, yrbc_int, zlbctype, zrbctype, &
                        zlbc_int, zrbc_int
                        
    use misc_module, only: to_upper
    implicit none
    integer, intent(in) :: ilo,ihi,jlo,jhi,klo,khi,nx,ny,nz

    ! is the block at the boundary
    at_xlboundary = .false.
    at_xrboundary = .false.
    at_ylboundary = .false.
    at_yrboundary = .false.
    at_zlboundary = .false.
    at_zrboundary = .false.

    if (ilo == 1)  at_xlboundary = .true.
    if (ihi == nx) at_xrboundary = .true.
#if NDIM >= 2
    if (jlo == 1)  at_ylboundary = .true.
    if (jhi == ny) at_yrboundary = .true.
#endif
#if NDIM == 3
    if (klo == 1)  at_zlboundary = .true.
    if (khi == nz) at_zrboundary = .true.
#endif

    ! defaults for inactive axes
    ylbc_int = FLOW
    yrbc_int = FLOW
    zlbc_int = FLOW
    zrbc_int = FLOW

    if (to_upper(trim(xlbctype)) == "FLOW") then 
        xlbc_int = FLOW
    else if (to_upper(trim(xlbctype)) == "REFLECT") then 
        xlbc_int = REFLECT
    else if (to_upper(trim(xlbctype)) == "PERIODIC") then
        xlbc_int = PERIODIC
    else if (to_upper(trim(xlbctype)) == "CUSTOM") then
        xlbc_int = CUSTOM
    else
        print *, "X-left BC options are: flow, reflect, periodic, custom "
        stop
    end if

     if (to_upper(trim(xrbctype)) == "FLOW") then
        xrbc_int = FLOW
    else if (to_upper(trim(xrbctype)) == "REFLECT") then 
        xrbc_int = REFLECT
    else if (to_upper(trim(xrbctype)) == "PERIODIC") then
        xrbc_int = PERIODIC
    else if (to_upper(trim(xrbctype)) == "CUSTOM") then
        xrbc_int = CUSTOM
    else
        print *, "X-right BC options are: flow, reflect, periodic, custom "
        stop
    end if

#if NDIM >= 2    
    if (to_upper(trim(ylbctype)) == "FLOW") then
        ylbc_int = FLOW
    else if (to_upper(trim(ylbctype)) == "REFLECT") then
        ylbc_int = REFLECT
    else if (to_upper(trim(ylbctype)) == "PERIODIC") then
        ylbc_int = PERIODIC
    else if (to_upper(trim(ylbctype)) == "CUSTOM") then
        ylbc_int = CUSTOM
    else
        print *, "Y-left BC options are: flow, reflect, periodic, custom "
        stop
    end if

     if (to_upper(trim(yrbctype)) == "FLOW") then
        yrbc_int = FLOW
    else if (to_upper(trim(yrbctype)) == "REFLECT") then 
        yrbc_int = REFLECT
    else if (to_upper(trim(yrbctype)) == "PERIODIC") then
        yrbc_int = PERIODIC
    else if (to_upper(trim(yrbctype)) == "CUSTOM") then
        yrbc_int = CUSTOM
    else
        print *, "Y-right BC options are: flow, reflect, periodic, custom "
        stop
    end if
#endif

#if NDIM == 3
    if (to_upper(trim(zlbctype)) == "FLOW") then
        zlbc_int = FLOW
    else if (to_upper(trim(zlbctype)) == "REFLECT") then
        zlbc_int = REFLECT
    else if (to_upper(trim(zlbctype)) == "PERIODIC") then
        zlbc_int = PERIODIC
    else if (to_upper(trim(zlbctype)) == "CUSTOM") then
        zlbc_int = CUSTOM
    else
        print *, "Z-left BC options are: flow, reflect, periodic, custom "
        stop
    end if

     if (to_upper(trim(zrbctype)) == "FLOW") then
        zrbc_int = FLOW
    else if (to_upper(trim(zrbctype)) == "REFLECT") then 
        zrbc_int = REFLECT
    else if (to_upper(trim(zrbctype)) == "PERIODIC") then
        zrbc_int = PERIODIC
    else if (to_upper(trim(zrbctype)) == "CUSTOM") then
        zrbc_int = CUSTOM
    else
        print *, "Z-right BC options are: flow, reflect, periodic, custom "
        stop
    end if
#endif

    !print *, "xrank, yrank, xl, xr, yl, yr = ", xrank, yrank, at_xlboundary, at_xrboundary, &
    !                                                          at_ylboundary, at_yrboundary
 end subroutine boundary_init

end module boundary_func
