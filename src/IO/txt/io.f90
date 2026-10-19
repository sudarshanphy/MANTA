module io_module
#include "param.h"

  implicit none
contains
  
  subroutine write_output(t, step, outputno, restart_init)
   use sim_data, only: nx, ny, nz, gamma, basenm, dx, dy, dz, &
                       mainVar, ilo, ihi, jlo, jhi, klo, khi, &
                       myrank, xblk, yblk, zblk, lnx, lny, lnz, &
                       iGlo, jGlo, kGlo, iGhi, jGhi, kGhi, outdir

   use grid_func, only: get_coords
   implicit none
   real, intent(in) :: t
   integer, intent(in) :: step, outputno
   real, pointer :: solnVar(:,:,:,:)
   logical, optional, intent(in) :: restart_init
   integer :: ionum, i, j, k
   character(len=256) :: fname
   character(len=10) :: int_to_str, rank_to_str

   real, allocatable:: x(:), y(:), z(:)

  allocate(x(iGlo:iGhi))
  allocate(y(jGlo:jGhi))
  allocate(z(kGlo:kGhi))

  call get_coords(IAXIS,iGlo,iGhi,x)
  call get_coords(JAXIS,jGlo,jGhi,y)
  call get_coords(KAXIS,kGlo,kGhi,z)

  solnVar(1:,iGlo:,jGlo:,kGlo:) => mainVar
  ionum = 10
  
12   format (1x, 50(es25.18, :, 1x))
  
  write(int_to_str, "(I4.4)") outputno
  write(rank_to_str, "(I4.4)") myrank

  if (present(restart_init)) then
    fname = trim(adjustl(outdir))//"/"//trim(adjustl(basenm))//"_reinit_"//trim(adjustl(rank_to_str))//&
      "_"//trim(adjustl(int_to_str))//".dat"
  else
    fname = trim(adjustl(outdir))//"/"//trim(adjustl(basenm))//"_"//trim(adjustl(rank_to_str))//&
      "_"//trim(adjustl(int_to_str))//".dat"
  endif

  open(unit=ionum, file=fname, status="replace")
  write(ionum, *) "#time = ", t, ", step = ", step
  write(ionum, *) "#nx = ", nx , ", dx = ", dx
  write(ionum, *) "#ny = ", ny , ", dy = ", dy
  write(ionum, *) "#lnx = ", lnx, ", xblk = ", xblk
  write(ionum, *) "#lny = ", lny, ", yblk = ", yblk
#if NDIM == 3
  write(ionum, *) "#nz = ", nz, ", dz = ", dz
  write(ionum, *) "#lnz = ", lnz, ", zblk = ", zblk
#endif
  write(ionum, *) "#gamma = ", gamma
#if NDIM != 2
  write(ionum, *) "#ndim = ", NDIM
#endif
#if NDIM == 1
  write(ionum, '(A)', advance='no') "# xcenter "
#elif NDIM == 2
  write(ionum, '(A)', advance='no') "# xcenter ycenter "
#else
  write(ionum, '(A)', advance='no') "# xcenter ycenter zcenter "
#endif
#ifdef MHD
  write(ionum, '(A)') "dens velx vely velz pres bmfx bmfy bmfz bpsi ener"
#else
  write(ionum, '(A)') "dens velx vely velz pres ener"
#endif
  ! keep the original x/y text ordering
  do k = klo, khi
    do i = ilo, ihi
      do j = jlo, jhi
#if NDIM == 1
        write(ionum, 12) x(i), solnVar(:,i,j,k)
#elif NDIM == 2
        write(ionum, 12) x(i), y(j), solnVar(:,i,j,k)
#else
        write(ionum, 12) x(i), y(j), z(k), solnVar(:,i,j,k)
#endif
      end do
      write(ionum, *) "####"
    end do
  end do
  close(ionum)
  deallocate(x,y,z)
  nullify(solnVar)
  end subroutine write_output
end module io_module
