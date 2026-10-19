module io_module
  use iso_c_binding, only: c_char, c_int, c_null_char
#include "param.h"
  implicit none
  interface
    function c_rename(old_name, new_name) bind(C, name="rename") result(status)
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: old_name(*), new_name(*)
      integer(c_int) :: status
    end function c_rename
  end interface
contains

  subroutine write_output(t, step, outputno, restart_init)
    use mpi, only: MPI_COMM_WORLD, MPI_Barrier, MPI_Comm_size
    use sim_data, only: nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, basenm, outdir, myrank, &
                        xblk, yblk, lnx, lny, lnz, &
                        ilo, ihi, jlo, jhi, klo, khi, iGlo, jGlo, kGlo, mainVar
    implicit none

    real,    intent(in) :: t
    integer, intent(in) :: step, outputno
    logical, optional, intent(in) :: restart_init

    real, pointer :: solnVar(:,:,:,:)
#ifdef MHD
    real(8), allocatable :: divB(:,:,:), divB_error(:,:,:)
#endif

    integer :: ierr, nprocs
    integer :: i0p, i1p, j0p, j1p, k0p, k1p
    integer :: whole_i0, whole_i1, whole_j0, whole_j1, whole_k0, whole_k1
    character(len=64) :: value_str
    integer :: ios, status
    logical :: is_restart
    integer :: ionum
    integer :: r, xr, yr, zr, ilo_r, ihi_r, jlo_r, jhi_r, klo_r, khi_r
    character(len=10)  :: int_to_str, rank_to_str
    character(len=:), allocatable :: subdir, vti_name, pvti_name, cmd, source_name
    real(8) :: yorigin, zorigin, dz_out, dy_out

    ! indices from param.h
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

    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    write(int_to_str,  "(I0.4)") outputno
    write(rank_to_str, "(I0.4)") myrank

    is_restart = .false.
    if (present(restart_init)) is_restart = restart_init

    if (is_restart) then
      subdir    = trim(adjustl(outdir))//"/vti_reinit_"//trim(adjustl(int_to_str))
      pvti_name = trim(adjustl(outdir))//"/"//trim(adjustl(basenm))//"_reinit_"//&
        trim(adjustl(int_to_str))//".pvti"
      vti_name  = trim(subdir)//"/"//trim(adjustl(basenm))//"_reinit_"//&
        trim(adjustl(rank_to_str))//"_"//trim(adjustl(int_to_str))//".vti"
    else
      subdir    = trim(adjustl(outdir))//"/vti_"//trim(adjustl(int_to_str))
      pvti_name = trim(adjustl(outdir))//"/"//trim(adjustl(basenm))//"_"//&
        trim(adjustl(int_to_str))//".pvti"
      vti_name  = trim(subdir)//"/"//trim(adjustl(basenm))//"_"//&
        trim(adjustl(rank_to_str))//"_"//trim(adjustl(int_to_str))//".vti"
    endif

    if (myrank == 0) then
      cmd = 'mkdir -p -- '//shell_quote(subdir)
      call execute_command_line(cmd, exitstat=status, cmdstat=ios)
      call check_io(ios, 'create output directory')
      call check_io(status, 'create output directory')
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! local point extents (inclusive)
    i0p = ilo - 1
    i1p = ihi
#if NDIM >= 2
    j0p = jlo - 1
    j1p = jhi
#else
    j0p = 0
    j1p = 0
#endif
#if NDIM == 3
    k0p = klo - 1
    k1p = khi
#else
    k0p = 0
    k1p = 0
#endif

    ! global point extents
    whole_i0 = 0
    whole_i1 = nx
#if NDIM >= 2
    whole_j0 = 0
    whole_j1 = ny
#else
    whole_j0 = 0
    whole_j1 = 0
#endif
    whole_k0 = 0
#if NDIM == 3
    whole_k1 = nz
#else
    whole_k1 = 0
#endif

    ! collapsed axes have zero extent
#if NDIM == 3
    zorigin = zmin
    dz_out = dz
#else
    zorigin = 0.0
    dz_out = 1.0
#endif
#if NDIM >= 2
    yorigin = ymin
    dy_out  = dy
#else
    yorigin = 0.0
    dy_out  = 1.0
#endif

    solnVar(1:, iGlo:, jGlo:, kGlo:) => mainVar
#ifdef MHD
    allocate(divB(ilo:ihi,jlo:jhi,klo:khi), divB_error(ilo:ihi,jlo:jhi,klo:khi))
    call get_divergence(divB, divB_error)
#endif

    ! ----------------------------
    ! Write per-rank .vti
    ! ----------------------------
    ionum = 10
    open(unit=ionum, file=vti_name//'.tmp', status="replace", action="write", &
         form="formatted", iostat=ios)
    call check_io(ios, 'open '//vti_name)

    call write_xml(ionum, '<?xml version="1.0"?>')
    call write_xml(ionum, '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">')
    ! each file contains only this rank's cells
    call write_xml(ionum, '  <ImageData WholeExtent="'//trim(extent_str(i0p,i1p,j0p,j1p,k0p,k1p))// &
                       '" Origin="'//trim(real3_str(xmin,yorigin,zorigin))// &
                       '" Spacing="'//trim(real3_str(dx,dy_out,dz_out))//'">')

    call write_xml(ionum, '      <FieldData>')
    call write_xml(ionum, '        <DataArray type="Int32" Name="CYCLE" NumberOfTuples="1" format="ascii">')
    write(value_str,'(I0)') step
    call write_xml(ionum, trim(value_str))
    call write_xml(ionum, '        </DataArray>')
    call write_xml(ionum, '        <DataArray type="Float64" Name="TIME" NumberOfTuples="1" format="ascii">')
    write(value_str,'(ES25.18)') t
    call write_xml(ionum, trim(value_str))
    call write_xml(ionum, '        </DataArray>')
    call write_xml(ionum, '      </FieldData>')
    
    call write_xml(ionum, '    <Piece Extent="'//trim(extent_str(i0p,i1p,j0p,j1p,k0p,k1p))//'">')

    call write_xml(ionum, '      <PointData/>')
    call write_xml(ionum, '      <CellData Scalars="dens">')

    call write_scalar(u=ionum, name="dens",  solnVar=solnVar, ivar=IDENS, ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)

    call write_scalar(u=ionum, name="velx",  solnVar=solnVar, ivar=IVX,   ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)
    call write_scalar(u=ionum, name="vely",  solnVar=solnVar, ivar=IVY,   ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)
    call write_scalar(u=ionum, name="velz",  solnVar=solnVar, ivar=IVZ,   ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)

    call write_scalar(u=ionum, name="pres",  solnVar=solnVar, ivar=IPRES, ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)

#ifdef MHD
    call write_scalar(u=ionum, name="bmfx",  solnVar=solnVar, ivar=IBX,   ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)
    call write_scalar(u=ionum, name="bmfy",  solnVar=solnVar, ivar=IBY,   ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)
    call write_scalar(u=ionum, name="bmfz",  solnVar=solnVar, ivar=IBZ,   ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)
    call write_scalar(u=ionum, name="bpsi",  solnVar=solnVar, ivar=IBPSI, ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)
#endif

    call write_scalar(u=ionum, name="ener",  solnVar=solnVar, ivar=IENER, ilo=ilo, ihi=ihi, &
                      jlo=jlo, jhi=jhi, klo=klo, khi=khi)

#ifdef MHD
    ! diagnostics follow the primitive fields
    call write_diagnostic(ionum, 'divB', divB)
    call write_diagnostic(ionum, 'divB_error', divB_error)
#endif

    call write_xml(ionum, '      </CellData>')
    call write_xml(ionum, '    </Piece>')
    call write_xml(ionum, '  </ImageData>')
    call write_xml(ionum, '</VTKFile>')

    close(ionum, iostat=ios)
    call check_io(ios, 'close '//vti_name)
    call publish_file(vti_name//'.tmp', vti_name)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! ----------------------------
    ! Write .pvti on rank 0 (deterministic extents)
    ! ----------------------------
    if (myrank == 0) then
      ionum = 20
      open(unit=ionum, file=pvti_name//'.tmp', status="replace", action="write", &
           form="formatted", iostat=ios)
      call check_io(ios, 'open '//pvti_name)
      call write_xml(ionum, '<?xml version="1.0"?>')
      call write_xml(ionum, '<VTKFile type="PImageData" version="0.1" byte_order="LittleEndian">')
      call write_xml(ionum, '  <PImageData WholeExtent="'// &
                     trim(extent_str(whole_i0,whole_i1,whole_j0,whole_j1,whole_k0,whole_k1))// &
                         '" GhostLevel="0" Origin="'//trim(real3_str(xmin,yorigin,zorigin))// &
                         '" Spacing="'//trim(real3_str(dx,dy_out,dz_out))//'">')

      !write(ionum,'(A)') '      <FieldData>'
      !write(ionum,'(A)') '        <DataArray type="Int32" Name="CYCLE" NumberOfTuples="1" format="ascii">'
      !write(ionum,'(I0)') step
      !write(ionum,'(A)') '        </DataArray>'
      !write(ionum,'(A)') '        <DataArray type="Float64" Name="TIME" NumberOfTuples="1" format="ascii">'
      !write(ionum,'(ES25.18)') dble(t)
      !write(ionum,'(A)') '        </DataArray>'
      !write(ionum,'(A)') '      </FieldData>'

      call write_xml(ionum, '    <PPointData/>')
      call write_xml(ionum, '    <PCellData Scalars="dens">')
      call write_p_scalar(ionum, "dens")

      call write_p_scalar(ionum, "velx")
      call write_p_scalar(ionum, "vely")
      call write_p_scalar(ionum, "velz")

      call write_p_scalar(ionum, "pres")

#ifdef MHD
      call write_p_scalar(ionum, "bmfx")
      call write_p_scalar(ionum, "bmfy")
      call write_p_scalar(ionum, "bmfz")
      call write_p_scalar(ionum, "bpsi")
#endif

      call write_p_scalar(ionum, "ener")
#ifdef MHD
      call write_p_scalar(ionum, "divB")
      call write_p_scalar(ionum, "divB_error")
#endif
      call write_xml(ionum, '    </PCellData>')

      do r = 0, nprocs-1
        xr = mod(r, xblk)
        yr = mod(r / xblk, yblk)
        zr = r / (xblk * yblk)

        ilo_r = xr * lnx + 1
        ihi_r = (xr + 1) * lnx
#if NDIM >= 2
        jlo_r = yr * lny + 1
        jhi_r = (yr + 1) * lny
#else
        jlo_r = 1
        jhi_r = 0
#endif

#if NDIM == 3
        klo_r = zr * lnz + 1
        khi_r = (zr + 1) * lnz
#else
        klo_r = 1
        khi_r = 0
#endif

        write(rank_to_str, "(I0.4)") r

        if (is_restart) then
          source_name = 'vti_reinit_'//trim(adjustl(int_to_str))//'/'//trim(adjustl(basenm))// &
                        '_reinit_'//trim(adjustl(rank_to_str))//'_'//trim(adjustl(int_to_str))//'.vti'
        else
          source_name = 'vti_'//trim(adjustl(int_to_str))//'/'//trim(adjustl(basenm))// &
                        '_'//trim(adjustl(rank_to_str))//'_'//trim(adjustl(int_to_str))//'.vti'
        endif
        call write_xml(ionum, '    <Piece Extent="'// &
                       trim(extent_str(ilo_r-1,ihi_r,jlo_r-1,jhi_r,klo_r-1,khi_r))// &
                       '" Source="'//xml_escape(source_name)//'"/>')
      end do

      call write_xml(ionum, '  </PImageData>')
      call write_xml(ionum, '</VTKFile>')

      close(ionum, iostat=ios)
      call check_io(ios, 'close '//pvti_name)
      call publish_file(pvti_name//'.tmp', pvti_name)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

#ifdef MHD
    deallocate(divB, divB_error)
#endif
    nullify(solnVar)

  end subroutine write_output


#ifdef MHD
  subroutine get_divergence(divB, divB_error)
    use sim_data, only: mainVar, ilo, ihi, jlo, jhi, klo, khi, dx, dy, dz, smallf
    real(8), intent(out) :: divB(ilo:ihi,jlo:jhi,klo:khi)
    real(8), intent(out) :: divB_error(ilo:ihi,jlo:jhi,klo:khi)
    real(8) :: widths(MAXDIM), h, Bmag
    integer :: i, j, k

    ! the caller has already filled guard cells
    widths = [dx,dy,dz]
    if (any(widths(1:NDIM) <= 0.0)) call check_io(1, 'nonpositive cell width for divB')
    h = minval(widths(1:NDIM))
    do k = klo, khi
      do j = jlo, jhi
        do i = ilo, ihi
          divB(i,j,k) = (mainVar(BMFX_VAR,i+1,j,k)-mainVar(BMFX_VAR,i-1,j,k))/(2.0*dx)
#if NDIM >= 2
          divB(i,j,k) = divB(i,j,k) + &
              (mainVar(BMFY_VAR,i,j+1,k)-mainVar(BMFY_VAR,i,j-1,k))/(2.0*dy)
#endif
#if NDIM == 3
          divB(i,j,k) = divB(i,j,k) + &
              (mainVar(BMFZ_VAR,i,j,k+1)-mainVar(BMFZ_VAR,i,j,k-1))/(2.0*dz)
#endif
          Bmag = sqrt(sum(mainVar(BMFX_VAR:BMFZ_VAR,i,j,k)**2))
          divB_error(i,j,k) = h*abs(divB(i,j,k))/max(Bmag,smallf)
        end do
      end do
    end do
  end subroutine get_divergence
#endif

  ! helper functions

  pure function extent_str(i0,i1,j0,j1,k0,k1) result(s)
    integer, intent(in) :: i0,i1,j0,j1,k0,k1
    character(len=128) :: s
    write(s,'(I0,1x,I0,1x,I0,1x,I0,1x,I0,1x,I0)') i0,i1,j0,j1,k0,k1
  end function extent_str

  pure function real3_str(a,b,c) result(s)
    real(8), intent(in) :: a,b,c
    character(len=128) :: s
    write(s,'(ES25.18,1x,ES25.18,1x,ES25.18)') a,b,c
  end function real3_str

  subroutine write_p_scalar(u, name)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    call write_xml(u, '      <PDataArray type="Float64" Name="'//trim(name)//'" NumberOfComponents="1"/>')
  end subroutine write_p_scalar

  subroutine write_scalar(u, name, solnVar, ivar, ilo, ihi, jlo, jhi, klo, khi)
    integer, intent(in) :: u, ivar, ilo, ihi, jlo, jhi, klo, khi
    character(len=*), intent(in) :: name
    real, pointer :: solnVar(:,:,:,:)
    integer :: i, j, k, ios

    call write_xml(u, '        <DataArray type="Float64" Name="'//trim(name)// &
                     '" NumberOfComponents="1" format="ascii">')
    ! VTK order: x first, then y, then z
    write(u,'(6(ES25.18,1X))', iostat=ios) &
        (((solnVar(ivar,i,j,k), i=ilo,ihi), j=jlo,jhi), k=klo,khi)
    call check_io(ios, 'write ASCII VTK field '//name)
    call write_xml(u, '        </DataArray>')
  end subroutine write_scalar

#ifdef MHD
  subroutine write_diagnostic(u, name, values)
    use sim_data, only: ilo, ihi, jlo, jhi, klo, khi
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    real(8), intent(in) :: values(ilo:ihi,jlo:jhi,klo:khi)
    integer :: ios

    call write_xml(u, '        <DataArray type="Float64" Name="'//trim(name)// &
                     '" NumberOfComponents="1" format="ascii">')
    write(u,'(6(ES25.18,1X))', iostat=ios) values
    call check_io(ios, 'write ASCII VTK diagnostic '//name)
    call write_xml(u, '        </DataArray>')
  end subroutine write_diagnostic
#endif

  subroutine write_xml(u, line)
    integer, intent(in) :: u
    character(len=*), intent(in) :: line
    character(len=11) :: form
    integer :: ios
    inquire(unit=u, form=form)
    if (form == 'UNFORMATTED') then
      write(u, iostat=ios) trim(line)//achar(10)
    else
      write(u,'(A)', iostat=ios) trim(line)
    end if
    call check_io(ios, 'write VTK header')
  end subroutine write_xml

  subroutine check_io(ios, message)
    use mpi, only: MPI_COMM_WORLD, MPI_Abort
    use sim_data, only: myrank
    integer, intent(in) :: ios
    character(len=*), intent(in) :: message
    integer :: ierr
    if (ios == 0) return
    print *, 'VTK I/O error: rank ', myrank, ': ', trim(message), ' status=', ios
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    error stop 'VTK I/O failed'
  end subroutine check_io

  subroutine publish_file(temporary, filename)
    character(len=*), intent(in) :: temporary, filename
    integer(c_int) :: status
    status = c_rename(trim(temporary)//c_null_char, trim(filename)//c_null_char)
    call check_io(status, 'publish '//filename)
  end subroutine publish_file

  function shell_quote(value) result(quoted)
    character(len=*), intent(in) :: value
    character(len=:), allocatable :: quoted
    integer :: i
    quoted = achar(39)
    do i = 1, len_trim(value)
      if (value(i:i) == achar(39)) then
        quoted = quoted//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
      else
        quoted = quoted//value(i:i)
      end if
    end do
    quoted = quoted//achar(39)
  end function shell_quote

  function xml_escape(value) result(escaped)
    character(len=*), intent(in) :: value
    character(len=:), allocatable :: escaped
    integer :: i
    escaped = ''
    do i = 1, len_trim(value)
      select case (value(i:i))
      case ('&')
        escaped = escaped//'&amp;'
      case ('"')
        escaped = escaped//'&quot;'
      case ('<')
        escaped = escaped//'&lt;'
      case ('>')
        escaped = escaped//'&gt;'
      case default
        escaped = escaped//value(i:i)
      end select
    end do
  end function xml_escape

end module io_module
