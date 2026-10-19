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
    character(len=256) :: fname, subdir
    integer :: u, ios
    integer :: ncells
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

    write(int_to_str,  "(I4.4)") fno
    write(rank_to_str, "(I4.4)") myrank

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
      if (ios /= 0) stop "restart_problem: cannot open .vti restart file"
    end if

    ncells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1)

    solnVar(1:, iGlo:, jGlo:, kGlo:) => mainVar
    solnVar = 0.0

    ! time
    call seek_to_dataarray(u, 'TIME')
    call read_n_reals_from_dataarray(u, 1, buf)
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
    integer, intent(in) :: u, ncells, ivar, ilo, ihi, jlo, jhi, klo, khi
    character(len=*), intent(in) :: name
    real, pointer :: solnVar(:,:,:,:)
    real(8), allocatable :: buf(:)
    integer :: i, j, k, idx
    character(len=512) :: header

    call seek_to_dataarray(u, name, header)
    if (index(header, 'format="appended"') > 0) then
      call read_binary_array(u, header, ncells, buf)
    else if (index(header, 'format="ascii"') > 0) then
      call read_n_reals_from_dataarray(u, ncells, buf)
    else
      error stop "restart_problem: unsupported VTK array format"
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
      if (ios /= 0) stop "seek_to_dataarray: DataArray not found"
      if (index(line, '</ImageData>') > 0) error stop "seek_to_dataarray: DataArray not found"
      if (index(line,'<DataArray') > 0 .and. index(line,'Name="'//trim(name)//'"') > 0) exit
    end do
    if (present(header)) header = line
  end subroutine seek_to_dataarray

  subroutine read_binary_array(u, header, n, arr)
    integer, intent(in) :: u, n
    character(len=*), intent(in) :: header
    real(8), allocatable, intent(out) :: arr(:)
    integer :: raw, ios, p, q, used
    integer(int64) :: offset, base, nbytes
    character(len=256) :: fname
    character(len=512) :: line
    character :: ch

    p = index(header, 'offset="')
    if (p == 0) error stop "restart_problem: missing binary offset"
    p = p + len('offset="')
    q = index(header(p:), '"')
    if (q <= 1) error stop "restart_problem: bad binary offset"
    read(header(p:p+q-2), *, iostat=ios) offset
    if (ios /= 0) error stop "restart_problem: bad binary offset"
    if (offset < 0) error stop "restart_problem: negative binary offset"

    inquire(unit=u, name=fname)
    close(u)
    open(newunit=raw, file=trim(fname), status="old", action="read", &
         access="stream", form="unformatted", iostat=ios)
    if (ios /= 0) error stop "restart_problem: cannot open binary stream"

    ! scan only the XML header, stopping before binary data
    line = ''
    used = 0
    do
      read(raw, iostat=ios) ch
      if (ios /= 0) error stop "restart_problem: missing appended data"
      if (ch == '>') then
        if (index(line, '<AppendedData') > 0) then
          if (index(line, 'encoding="raw"') == 0) error stop "restart_problem: expected raw encoding"
          exit
        end if
      end if
      if (ch == achar(10)) then
        if (index(line, '<VTKFile') > 0) then
          if (index(line, 'header_type="UInt64"') == 0) error stop "restart_problem: expected UInt64 header"
          if (index(line, 'byte_order="LittleEndian"') == 0) error stop "restart_problem: expected little endian"
          if (index(line, 'compressor=') > 0) error stop "restart_problem: compressed arrays unsupported"
        end if
        line = ''
        used = 0
      else
        used = used + 1
        if (used > len(line)) error stop "restart_problem: XML header line too long"
        line(used:used) = ch
      end if
    end do

    ! offsets start immediately after the underscore
    do
      read(raw, iostat=ios) ch
      if (ios /= 0) error stop "restart_problem: missing appended marker"
      if (ch > ' ') exit
    end do
    if (ch /= '_') error stop "restart_problem: bad appended marker"
    inquire(unit=raw, pos=base)
    read(raw, pos=base+offset, iostat=ios) nbytes
    if (ios /= 0) error stop "restart_problem: cannot read binary size"
    if (nbytes /= 8_int64*int(n,int64)) error stop "restart_problem: binary size mismatch"
    allocate(arr(n))
    read(raw, iostat=ios) arr
    if (ios /= 0) error stop "restart_problem: incomplete binary field"
    close(raw)
    open(unit=u, file=trim(fname), status="old", action="read", form="formatted", iostat=ios)
    if (ios /= 0) error stop "restart_problem: cannot reopen XML header"
  end subroutine read_binary_array

  subroutine read_n_reals_from_dataarray(u, n, arr)
    integer, intent(in) :: u, n
    real(8), allocatable, intent(out) :: arr(:)
    character(len=512) :: line
    real(8) :: vals(256)
    integer :: ios, got, m, ii

    allocate(arr(n))
    got = 0

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) stop "read_n_reals_from_dataarray: unexpected EOF"
      if (index(line,'</DataArray>') > 0) exit

      call parse_reals_from_line(line, vals, m)
      do ii = 1, m
        got = got + 1
        if (got <= n) arr(got) = vals(ii)
      end do

      if (got >= n) then
        do
          if (index(line,'</DataArray>') > 0) exit
          read(u,'(A)', iostat=ios) line
          if (ios /= 0) exit
        end do
        exit
      end if
    end do

    if (got /= n) stop "read_n_reals_from_dataarray: wrong number of values"
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
      do while (p <= L .and. line(p:p) <= ' ')
        p = p + 1
      end do
      if (p > L) exit

      q = p
      do while (q <= L .and. line(q:q) > ' ')
        q = q + 1
      end do

      tok = line(p:q-1)

      if (index(tok,'<') == 0 .and. index(tok,'>') == 0) then
        read(tok,*,iostat=ios) tmp
        if (ios == 0) then
          nout = nout + 1
          if (nout <= vmax) vals(nout) = tmp
        end if
      end if

      p = q + 1
      if (nout >= vmax) exit
    end do
  end subroutine parse_reals_from_line

end module io_restart
