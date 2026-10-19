module grid_func
#include "param.h"

  implicit none
  contains

  subroutine grid_init()
    use sim_data
    use misc_module, only: to_upper
    use boundary_func, only: boundary_init
    implicit none

! inactive axes
#if NDIM == 1
    ny = 1
    yblk = 1
#endif
#if NDIM < 3
    nz = 1
    zblk = 1
#endif

    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    dz = (zmax - zmin) / nz
   
    ! xrank and yrank 
    xrank = mod(myrank, xblk)
    yrank = mod(myrank/xblk, yblk)
    zrank = myrank / (xblk * yblk)
    
    ! default guard cells
    Gpts = 2  
    !if (to_upper(trim(recon_method)) == "WENO3") Gpts = 2
    if (to_upper(trim(recon_method)) == "WENO5") Gpts = 3
    if (to_upper(trim(recon_method)) == "PPM") Gpts = 3

    nghost = 0
    nghost(1:NDIM) = Gpts

    ! local interior points 
    lnx = nx / xblk
    lny = ny / yblk
    lnz = nz / zblk
  
    ! compute index of interior points 
    ilo = xrank * lnx + 1
    ihi = (xrank + 1) * lnx  
    jlo = yrank * lny + 1
    jhi = (yrank + 1) * lny 
    klo = zrank * lnz + 1
    khi = (zrank + 1) * lnz 
 
    ! compute index of interior + guard cells
    iGlo = ilo - nghost(IAXIS) 
    iGhi = ihi + nghost(IAXIS) 
    jGlo = jlo - nghost(JAXIS) 
    jGhi = jhi + nghost(JAXIS) 
    kGlo = klo - nghost(KAXIS) 
    kGhi = khi + nghost(KAXIS) 
 
    ! total number of local points 
    lxTpts = lnx + 2 * nghost(IAXIS) 
    lyTpts = lny + 2 * nghost(JAXIS) 
    lzTpts = lnz + 2 * nghost(KAXIS) 
  
    ! local min and max
    lxmin = xmin + (ilo - 1) * dx
    lxmax = xmin + ihi * dx
    lymin = ymin + (jlo - 1) * dy
    lymax = ymin + jhi * dy
    lzmin = zmin + (klo - 1) * dz
    lzmax = zmin + khi * dz
  
    allocate(mainVar(NVAR_NUMBER, iGlo:iGhi, jGlo:jGhi, kGlo:kGhi))

    !set all the small values
    smalld = 1.0e-10
    smallp = 1.0e-10
    smalle = 1.0e-10
    small = 1.0e-100
 
    !print *, "myrank, xl, xr, yl, yr = ", myrank, at_xlboundary, at_xrboundary, &
    !                                              at_ylboundary, at_yrboundary 
    !print *, "myrank, ilo, ihi, jlo, jhi = ", myrank, ilo, ihi, jlo, jhi

    call boundary_init(ilo,ihi,jlo,jhi,klo,khi,nx,ny,nz)

  end subroutine grid_init

  subroutine grid_finalize()
    use sim_data, only: mainVar
    implicit none

    deallocate(mainVar)
  end subroutine

  subroutine get_coords(Dir,lo,hi,array)
    use sim_data, only: xmin, ymin, zmin, dx, dy, dz
    implicit none
    integer, intent(in) :: Dir
    integer, intent(in) :: lo, hi
    real(8), dimension(lo:hi), intent(out) :: array
    integer :: ii
    real(8) :: delta, lmin
    
    array(lo:hi) = 0.0
 
    select case (Dir)
    case (IAXIS)
      lmin = xmin
      delta = dx
    case (JAXIS)
      lmin = ymin
      delta = dy
    case (KAXIS)
      lmin = zmin
      delta = dz
    case default
      print *, "Direction is X, Y or Z"
      stop
    end select

    do ii = lo, hi
      array(ii) = lmin + (2 * ii - 1) * delta/2.0
    end do 

  end subroutine  

end module grid_func
