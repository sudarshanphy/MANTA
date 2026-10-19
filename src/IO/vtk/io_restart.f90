module io_restart
  use iso_fortran_env, only: int64
#include "param.h"
  implicit none
contains

  subroutine restart_problem(fno, time)
    use sim_data, only: basenm, outdir, myrank, &
                        ilo, ihi, jlo, jhi, klo, khi, iGlo, jGlo, kGlo, mainVar
    implicit none

    integer, intent(in) :: fno
    real,    intent(out):: time

    real, pointer :: solnVar(:,:,:,:)

    character(len=10)  :: int_to_str, rank_to_str
    character(len=:), allocatable :: fname, subdir
    integer :: u, ios
    integer(int64) :: ncells
    real(8), allocatable :: buf(:)

    integer, parameter :: IDENS = DENS_VAR
    integer, parameter :: IVX   = VELX_VAR
    integer, parameter :: IVY   = VELY_VAR
    integer, parameter :: IVZ   = VELZ_VAR
    integer, parameter :: IPRES = PRES_VAR
    integer, parameter :: IENER = ENER_VAR
#ifdef MHD
    integer, parameter :: IBX   = BMFX_VAR
    integer, parameter :: IBY   = BMFY_VAR
    integer, parameter :: IBZ   = BMFZ_VAR
    integer, parameter :: IBPSI = BPSI_VAR
#endif

    write(int_to_str,  "(I0.4)") fno
    write(rank_to_str, "(I0.4)") myrank

    subdir = trim(adjustl(outdir))//"/vti_"//trim(adjustl(int_to_str))
    fname  = trim(subdir)//"/"//trim(adjustl(basenm))//"_"//trim(adjustl(rank_to_str))//"_"// &
             trim(adjustl(int_to_str))//".vti"

    u = 56
    open(unit=u, file=fname, status="old", action="read", form="formatted", iostat=ios)
    if (ios /= 0) then
      subdir = trim(adjustl(outdir))//"/vti_reinit_"//trim(adjustl(int_to_str))
      fname  = trim(subdir)//"/"//trim(adjustl(basenm))//"_reinit_"//trim(adjustl(rank_to_str))//"_"// &
               trim(adjustl(int_to_str))//".vti"
      open(unit=u, file=fname, status="old", action="read", form="formatted", iostat=ios)
      if (ios /= 0) call io_fail('restart_problem: cannot open .vti restart file')
    end if

    call validate_restart(u)
    ncells = int(ihi-ilo+1,int64) * int(jhi-jlo+1,int64) * int(khi-klo+1,int64)

    solnVar(1:, iGlo:, jGlo:, kGlo:) => mainVar
    solnVar = 0.0

    ! time
    call seek_to_dataarray(u, 'TIME')
    call read_n_reals_from_dataarray(u, 1_int64, buf)
    time = real(buf(1))
    deallocate(buf)

    call read_scalar_into(u, "dens",  ncells, solnVar, IDENS, ilo, ihi, jlo, jhi, klo, khi)

    call read_scalar_into(u, "velx",  ncells, solnVar, IVX,   ilo, ihi, jlo, jhi, klo, khi)
    call read_scalar_into(u, "vely",  ncells, solnVar, IVY,   ilo, ihi, jlo, jhi, klo, khi)
    call read_scalar_into(u, "velz",  ncells, solnVar, IVZ,   ilo, ihi, jlo, jhi, klo, khi)

    call read_scalar_into(u, "pres",  ncells, solnVar, IPRES, ilo, ihi, jlo, jhi, klo, khi)

#ifdef MHD
    call read_scalar_into(u, "bmfx",  ncells, solnVar, IBX,   ilo, ihi, jlo, jhi, klo, khi)
    call read_scalar_into(u, "bmfy",  ncells, solnVar, IBY,   ilo, ihi, jlo, jhi, klo, khi)
    call read_scalar_into(u, "bmfz",  ncells, solnVar, IBZ,   ilo, ihi, jlo, jhi, klo, khi)
    call read_scalar_into(u, "bpsi",  ncells, solnVar, IBPSI, ilo, ihi, jlo, jhi, klo, khi)
#endif

    call read_scalar_into(u, "ener",  ncells, solnVar, IENER, ilo, ihi, jlo, jhi, klo, khi)

    close(u)
    nullify(solnVar)

  end subroutine restart_problem


  ! -----------------------------
  ! Read a scalar DataArray into solnVar(ivar,:,:,:)
  ! -----------------------------
  subroutine read_scalar_into(u, name, ncells, solnVar, ivar, ilo, ihi, jlo, jhi, klo, khi)
    integer, intent(in) :: u, ivar, ilo, ihi, jlo, jhi, klo, khi
    integer(int64), intent(in) :: ncells
    character(len=*), intent(in) :: name
    real, pointer :: solnVar(:,:,:,:)
    real(8), allocatable :: buf(:)
    integer :: i, j, k
    integer(int64) :: idx
    character(len=512) :: header

    call seek_to_dataarray(u, name, header)
    if (index(header, 'type="Float64"') == 0) call io_fail('expected Float64 field '//name)
    if (index(header, 'format="appended"') > 0) then
      call read_binary_array(u, header, ncells, buf)
    else if (index(header, 'format="ascii"') > 0) then
      call read_n_reals_from_dataarray(u, ncells, buf)
    else
      call io_fail('restart_problem: unsupported VTK array format')
    end if

    idx = 0
    do k = klo, khi
      do j = jlo, jhi
        do i = ilo, ihi
          idx = idx + 1
          solnVar(ivar,i,j,k) = real(buf(idx))
        end do
      end do
    end do
    deallocate(buf)
  end subroutine read_scalar_into


  ! -----------------------------
  ! XML parsing helpers
  ! -----------------------------
  subroutine seek_to_dataarray(u, name, header)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(out) :: header
    character(len=512) :: line
    integer :: ios
    rewind(u)
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) call io_fail('seek_to_dataarray: DataArray not found')
      if (index(line, '</ImageData>') > 0) call io_fail('seek_to_dataarray: DataArray not found')
      if (index(line,'<DataArray') > 0 .and. index(line,'Name="'//trim(name)//'"') > 0) exit
    end do
    if (present(header)) header = line
  end subroutine seek_to_dataarray

  subroutine read_binary_array(u, header, n, arr)
    integer, intent(in) :: u
    integer(int64), intent(in) :: n
    character(len=*), intent(in) :: header
    real(8), allocatable, intent(out) :: arr(:)
    integer :: raw, ios, p, q, used
    integer(int64) :: offset, base, nbytes, filesize
    character(len=1024) :: fname
    character(len=512) :: line
    character :: ch

    p = index(header, 'offset="')
    if (p == 0) call io_fail('restart_problem: missing binary offset')
    p = p + len('offset="')
    q = index(header(p:), '"')
    if (q <= 1) call io_fail('restart_problem: bad binary offset')
    read(header(p:p+q-2), *, iostat=ios) offset
    if (ios /= 0) call io_fail('restart_problem: bad binary offset')
    if (offset < 0) call io_fail('restart_problem: negative binary offset')

    inquire(unit=u, name=fname)
    close(u)
    open(newunit=raw, file=trim(fname), status="old", action="read", &
         access="stream", form="unformatted", iostat=ios)
    if (ios /= 0) call io_fail('restart_problem: cannot open binary stream')

    ! scan only the XML header, stopping before binary data
    line = ''
    used = 0
    do
      read(raw, iostat=ios) ch
      if (ios /= 0) call io_fail('restart_problem: missing appended data')
      if (ch == '>') then
        if (index(line, '<AppendedData') > 0) then
          if (index(line, 'encoding="raw"') == 0) call io_fail('restart_problem: expected raw encoding')
          exit
        end if
      end if
      if (ch == achar(10)) then
        if (index(line, '<VTKFile') > 0) then
          if (index(line, 'header_type="UInt64"') == 0) call io_fail('restart_problem: expected UInt64 header')
          if (index(line, 'byte_order="LittleEndian"') == 0) call io_fail('restart_problem: expected little endian')
          if (index(line, 'compressor=') > 0) call io_fail('restart_problem: compressed arrays unsupported')
        end if
        line = ''
        used = 0
      else
        used = used + 1
        if (used > len(line)) call io_fail('restart_problem: XML header line too long')
        line(used:used) = ch
      end if
    end do

    ! offsets start immediately after the underscore
    do
      read(raw, iostat=ios) ch
      if (ios /= 0) call io_fail('restart_problem: missing appended marker')
      if (ch > ' ') exit
    end do
    if (ch /= '_') call io_fail('restart_problem: bad appended marker')
    inquire(unit=raw, pos=base)
    inquire(unit=raw, size=filesize)
    if (offset > filesize-base-7_int64) call io_fail('binary offset beyond end of file')
    read(raw, pos=base+offset, iostat=ios) nbytes
    if (ios /= 0) call io_fail('restart_problem: cannot read binary size')
    if (nbytes /= 8_int64*int(n,int64)) call io_fail('restart_problem: binary size mismatch')
    if (nbytes > filesize-base-offset-7_int64) call io_fail('incomplete binary field')
    allocate(arr(n))
    read(raw, iostat=ios) arr
    if (ios /= 0) call io_fail('restart_problem: incomplete binary field')
    close(raw)
    open(unit=u, file=trim(fname), status="old", action="read", form="formatted", iostat=ios)
    if (ios /= 0) call io_fail('restart_problem: cannot reopen XML header')
  end subroutine read_binary_array

  subroutine read_n_reals_from_dataarray(u, n, arr)
    integer, intent(in) :: u
    integer(int64), intent(in) :: n
    real(8), allocatable, intent(out) :: arr(:)
    character(len=512) :: line
    real(8) :: vals(256)
    integer :: ios, m, ii
    integer(int64) :: got

    allocate(arr(n))
    got = 0

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) call io_fail('read_n_reals_from_dataarray: unexpected EOF')
      if (index(line,'</DataArray>') > 0) exit

      call parse_reals_from_line(line, vals, m)
      if (got + int(m,int64) > n) call io_fail('too many values in ASCII field')
      do ii = 1, m
        got = got + 1
        arr(got) = vals(ii)
      end do
    end do

    if (got /= n) call io_fail('read_n_reals_from_dataarray: wrong number of values')
  end subroutine read_n_reals_from_dataarray

  subroutine parse_reals_from_line(line, vals, nout)
    character(len=*), intent(in) :: line
    real(8), intent(out) :: vals(:)
    integer, intent(out) :: nout

    integer :: L, p, q, ios, vmax
    character(len=64) :: tok
    real(8) :: tmp

    vmax = size(vals)
    vals = 0.0
    nout = 0

    L = len_trim(line)
    p = 1

    do while (p <= L)
      do while (p <= L)
        if (line(p:p) > ' ') exit
        p = p + 1
      end do
      if (p > L) exit

      q = p
      do while (q <= L)
        if (line(q:q) <= ' ') exit
        q = q + 1
      end do

      if (q-p > len(tok)) call io_fail('ASCII token too long')
      tok = line(p:q-1)

      if (index(tok,'<') == 0 .and. index(tok,'>') == 0) then
        read(tok,*,iostat=ios) tmp
        if (ios /= 0) call io_fail('invalid number in ASCII field')
        if (nout == vmax) call io_fail('too many numbers on ASCII line')
        nout = nout + 1
        vals(nout) = tmp
      end if

      p = q + 1

    end do
  end subroutine parse_reals_from_line

  subroutine validate_restart(u)
    use sim_data, only: ilo, ihi, jlo, jhi, klo, khi
    integer, intent(in) :: u
    integer :: expected(6), actual(6), p, q, ios, raw, count
    integer(int64) :: filesize
    character(len=1024) :: fname
    character(len=512) :: line
    character(len=:), allocatable :: tail

    ! reject unfinished files before reading fields
    inquire(unit=u, name=fname)
    close(u)
    open(newunit=raw, file=trim(fname), status='old', action='read', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) call io_fail('cannot open restart stream')
    inquire(unit=raw, size=filesize)
    count = int(min(filesize,128_int64))
    if (count == 0) call io_fail('empty restart file')
    allocate(character(len=count) :: tail)
    read(raw, pos=filesize-count+1, iostat=ios) tail
    if (ios /= 0) call io_fail('cannot read restart footer')
    if (index(tail, '</VTKFile>') == 0) call io_fail('incomplete restart file')
    close(raw)
    open(unit=u, file=trim(fname), status='old', action='read', form='formatted', iostat=ios)
    if (ios /= 0) call io_fail('cannot open restart header')

    expected = [ilo-1,ihi,0,0,0,0]
#if NDIM >= 2
    expected(3:4) = [jlo-1,jhi]
#endif
#if NDIM == 3
    expected(5:6) = [klo-1,khi]
#endif
    ! the piece must match this rank's interior
    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) call io_fail('missing restart Piece')
      if (index(line, '<Piece ') > 0) exit
      if (index(line, '</ImageData>') > 0) call io_fail('missing restart Piece')
    end do
    p = index(line, 'Extent="')
    if (p == 0) call io_fail('missing restart extent')
    p = p + len('Extent="')
    q = index(line(p:), '"')
    if (q <= 1) call io_fail('bad restart extent')
    read(line(p:p+q-2), *, iostat=ios) actual
    if (ios /= 0) call io_fail('bad restart extent')
    if (any(actual /= expected)) call io_fail('restart extent does not match grid/decomposition')
    rewind(u)
  end subroutine validate_restart

  subroutine io_fail(message)
    use mpi, only: MPI_COMM_WORLD, MPI_Abort
    use sim_data, only: myrank
    character(len=*), intent(in) :: message
    integer :: ierr
    print *, 'VTK restart error: rank ', myrank, ': ', trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    error stop 'VTK restart failed'
  end subroutine io_fail

end module io_restart
