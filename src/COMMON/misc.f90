module misc_module
#include "param.h"
  use mpi
  implicit none
contains

    function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page
    
         implicit none
    
         character(len=*), intent(in) :: strIn
         character(len=len(strIn)) :: strOut
         integer :: i,j
    
         do i = 1, len(strIn)
              j = iachar(strIn(i:i))
              if (j>= iachar("a") .and. j<=iachar("z") ) then
                   strOut(i:i) = achar(iachar(strIn(i:i))-32)
              else
                   strOut(i:i) = strIn(i:i)
              end if
         end do
    
    end function to_upper

    subroutine get_dt(dt)
      use sim_data, only: gamma, dx, dy, dz, cfl, mainVar, &
                          iGlo, iGhi, jGlo, jGhi, kGlo, kGhi, &
                          comm, ierr
      use mpi
#ifdef MHD
      use sim_data, only: ch
#endif
      implicit none
      real(8), dimension(iGlo:iGhi,jGlo:jGhi,kGlo:kGhi) :: xsmax, ysmax, zsmax
      real(8) :: cs, cs2, xcmax, ycmax, zcmax   !sound speed
#ifdef MHD
      real(8) :: cax, cay, caz, cfx, cfy, cfz, B2, cB2 !alfven wave and magnetosonic wave
      real(8) :: localch
#endif
      real(8), intent(out) :: dt
      real(8) :: localdt, rate, cell_width(MAXDIM)
      integer :: i, j, k
      real(8) :: max_xsmax, max_ysmax, max_zsmax
      real(8), pointer :: solnVar(:,:,:,:)

      cell_width = [dx,dy,dz]
      if (any(cell_width(1:NDIM) <= 0.0) .or. cfl <= 0.0) then
        print *, "Active cell widths and CFL must be positive"
        call MPI_Abort(comm, 1, ierr)
      end if

      xsmax = 0.0
      ysmax = 0.0
      zsmax = 0.0
      solnVar(1:,iGlo:,jGlo:,kGlo:) => mainVar(:,:,:,:)

      do k = kGlo, kGhi
        do j = jGlo, jGhi
          do i = iGlo, iGhi
            if (solnVar(DENS_VAR,i,j,k) <= 0.0 .or. solnVar(PRES_VAR,i,j,k) < 0.0) then
              print *, "Invalid density or pressure at i,j,k = ", i,j,k
              call MPI_Abort(comm, 1, ierr)
            end if
            cs2 = gamma * solnVar(PRES_VAR,i,j,k) / solnVar(DENS_VAR,i,j,k)

            if (cs2 < 0.0) then
              print *, "Imaginary sound speed!"
              print *, "At i,j,k, with gamma, pres, dens = ", i, j, k, gamma, &
                            solnVar(PRES_VAR,i,j,k), solnVar(DENS_VAR,i,j,k)
              call MPI_Abort(comm, 1, ierr)
            end if
            cs = sqrt(cs2)
            xcmax = cs
#if NDIM >= 2
            ycmax = cs
#endif
#if NDIM == 3
            zcmax = cs
#endif
#ifdef MHD
            cax = solnVar(BMFX_VAR,i,j,k) / sqrt(solnVar(DENS_VAR,i,j,k))
            B2 = sum(solnVar(BMFX_VAR:BMFZ_VAR,i,j,k)*solnVar(BMFX_VAR:BMFZ_VAR,i,j,k))
            cB2 = B2 / solnVar(DENS_VAR,i,j,k)
            cfx = sqrt(0.5 * ((cs2 + cB2) + &
                  sqrt(max((cs2 + cB2)**2 - 4.0*cs2*cax*cax, 0.0))))
            xcmax = cfx
#if NDIM >= 2
            cay = solnVar(BMFY_VAR,i,j,k) / sqrt(solnVar(DENS_VAR,i,j,k))
            cfy = sqrt(0.5 * ((cs2 + cB2) + &
                  sqrt(max((cs2 + cB2)**2 - 4.0*cs2*cay*cay, 0.0))))
            ycmax = cfy
#endif
#if NDIM == 3
            caz = solnVar(BMFZ_VAR,i,j,k) / sqrt(solnVar(DENS_VAR,i,j,k))
            cfz = sqrt(0.5 * ((cs2 + cB2) + &
                  sqrt(max((cs2 + cB2)**2 - 4.0*cs2*caz*caz, 0.0))))
            zcmax = cfz
#endif
#endif
            xsmax(i,j,k) = xcmax + abs(solnVar(VELX_VAR,i,j,k))
#if NDIM >= 2
            ysmax(i,j,k) = ycmax + abs(solnVar(VELY_VAR,i,j,k))
#endif
#if NDIM == 3
            zsmax(i,j,k) = zcmax + abs(solnVar(VELZ_VAR,i,j,k))
#endif
          end do
        end do
      end do

      max_xsmax = maxval(xsmax)
      max_ysmax = maxval(ysmax)
      max_zsmax = maxval(zsmax)

#ifdef MHD
      !speed for divergence cleaning
      localch = max(max_xsmax, max_ysmax, max_zsmax)
      ! get max ch from all the cores
      call MPI_ALLREDUCE(localch, ch, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

      ! cleaning waves propagate in every active direction
      max_xsmax = max(max_xsmax,ch)
#if NDIM >= 2
      max_ysmax = max(max_ysmax,ch)
#endif
#if NDIM == 3
      max_zsmax = max(max_zsmax,ch)
#endif
#endif

      !print *, "max xsmax, ysmax = ", max_xsmax, max_ysmax
      ! sum active wave rates for the unsplit RK update
      rate = max_xsmax / dx
#if NDIM >= 2
      rate = rate + max_ysmax / dy
#endif
#if NDIM == 3
      rate = rate + max_zsmax / dz
#endif
      localdt = huge(localdt)
      if (rate > 0.0) localdt = cfl / rate

      !print *, "local dt = ", localdt
      ! get min dt from all the cores
      call MPI_ALLREDUCE(localdt, dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
      !print *, "global dt = ", dt
      nullify(solnVar)

    end subroutine get_dt

    function compute_maxdivB(Bx, By, Bz, normalized, b_floor) result(maxdivB)
      use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                          iGlo, jGlo, kGlo, dx, dy, dz, smallf, comm, ierr
      implicit none
      real(8), intent(in), dimension(iGlo:,jGlo:,kGlo:) :: Bx, By, Bz
      logical, intent(in), optional :: normalized
      real(8), intent(in), optional :: b_floor
      real(8) :: maxdivB, localmax, divB, error, Bmag, floorB, h
      real(8) :: cell_width(MAXDIM)
      logical :: use_normalized
      integer :: i, j, k

      ! call on all ranks after guard cells have been filled
      use_normalized = .false.
      if (present(normalized)) use_normalized = normalized
      floorB = smallf
      if (present(b_floor)) floorB = b_floor
      if (use_normalized .and. floorB <= 0.0) then
        print *, "Magnetic normalization floor must be positive"
        call MPI_Abort(comm, 1, ierr)
      end if
      cell_width = [dx,dy,dz]
      if (any(cell_width(1:NDIM) <= 0.0)) then
        print *, "Active cell widths must be positive"
        call MPI_Abort(comm, 1, ierr)
      end if
      h = minval(cell_width(1:NDIM))

      localmax = 0.0
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            divB = (Bx(i+1,j,k) - Bx(i-1,j,k)) / (2.0 * dx)
#if NDIM >= 2
            divB = divB + (By(i,j+1,k) - By(i,j-1,k)) / (2.0 * dy)
#endif
#if NDIM == 3
            divB = divB + (Bz(i,j,k+1) - Bz(i,j,k-1)) / (2.0 * dz)
#endif
            !if (divB > tol) then
            !  print *, "High error in divergence of B (> 1.0e-4)!"
            !  !print *, "Error at x, y = ", xvals(i), yvals(j)
            !  print *, "Error = ", divB
            !  stop
            !end if
            error = abs(divB)
            if (use_normalized) then
              Bmag = sqrt(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)
              error = h * error / max(Bmag,floorB)
            end if
            localmax = max(localmax,error)
          end do
        end do
      end do

      call MPI_ALLREDUCE(localmax, maxdivB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    end function compute_maxdivB

end module misc_module
