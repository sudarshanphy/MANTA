module io_restart
#include "param.h"
  implicit none
contains
  subroutine restart_problem(fno, time)
    use sim_data, only: basenm, ilo, ihi, jlo, jhi, klo, khi, &
                        mainVar, myrank, outdir
    implicit none
    integer, intent(in) :: fno
    real, intent(out) :: time
    real, pointer :: solnVar(:,:,:,:)
    integer :: i, j, k, ionum, pos1, pos2, ios, file_ndim
    character(len=2048) :: line
    character(len=256) :: fname
    character(len=10) :: int_to_str, rank_to_str
    real :: coords(NDIM)
    logical :: have_time

    solnVar(1:,ilo:,jlo:,klo:) => mainVar(:,ilo:ihi,jlo:jhi,klo:khi)
    solnVar = 0.0

    ionum = 56
    write(int_to_str, "(I4.4)") fno
    write(rank_to_str, "(I4.4)") myrank

    fname = trim(adjustl(outdir))//"/"//trim(adjustl(basenm))//"_"//trim(adjustl(rank_to_str))//&
      "_"//trim(adjustl(int_to_str))//".dat"

    open(unit=ionum, file=fname, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      fname = trim(adjustl(outdir))//"/"//trim(adjustl(basenm))//"_reinit_"//trim(adjustl(rank_to_str))//&
        "_"//trim(adjustl(int_to_str))//".dat"
      open(unit=ionum, file=fname, status="old", action="read", iostat=ios)
      if (ios /= 0) error stop "restart_problem: cannot open text restart file"
    end if

    ! read headers up to the first cell
    have_time = .false.
    do
      call next_line()
      if (line(1:1) /= '#') exit
      if (index(line, '#time') == 1) then
        ! Find the first instance of whitespace.  Split label and data.
        pos1 = scan(line, '=')
        pos2 = scan(line, ',')
        if (pos1 == 0 .or. pos2 <= pos1) error stop "restart_problem: bad time header"
        read(line(pos1+1:pos2-1), *, iostat=ios) time
        if (ios /= 0) error stop "restart_problem: bad time value"
        have_time = .true.
      else if (index(line, '#ndim') == 1) then
        pos1 = scan(line, '=')
        read(line(pos1+1:), *, iostat=ios) file_ndim
        if (ios /= 0) error stop "restart_problem: bad dimension header"
        if (file_ndim /= NDIM) error stop "restart_problem: dimension mismatch"
      end if
    end do
    if (.not. have_time) error stop "restart_problem: missing time header"

    do k = klo, khi
      do i = ilo, ihi
        do j = jlo, jhi
          read(line, *, iostat=ios) coords, solnVar(NVAR_BEGIN:NVAR_END,i,j,k)
          if (ios /= 0) error stop "restart_problem: bad cell data"
          ! skip blank lines and row separators
          if (i == ihi .and. j == jhi .and. k == khi) exit
          do
            call next_line()
            if (line(1:1) /= '#') exit
          end do
        end do
      end do
    end do
    close(ionum)
    nullify(solnVar)

  contains
    subroutine next_line()
      do
        read(ionum, '(A)', iostat=ios) line
        if (ios /= 0) error stop "restart_problem: unexpected end of text file"
        line = adjustl(line)
        if (len_trim(line) > 0) exit
      end do
    end subroutine next_line
  end subroutine restart_problem
end module io_restart
