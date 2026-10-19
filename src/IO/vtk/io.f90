module io_module
  use iso_fortran_env, only: int64, real64
#include "param.h"
  implicit none
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

    integer :: ierr, nprocs
    integer :: i0p, i1p, j0p, j1p, k0p, k1p
    integer :: whole_i0, whole_i1, whole_j0, whole_j1, whole_k0, whole_k1
    integer :: ivar
    integer(int64) :: nbytes
    real(real64), allocatable :: field(:,:,:)
    character(len=64) :: value_str
    integer :: ionum
    integer :: r, xr, yr, zr, ilo_r, ihi_r, jlo_r, jhi_r, klo_r, khi_r
    character(len=10)  :: int_to_str, rank_to_str
    character(len=256) :: subdir, vti_name, pvti_name, cmd
    real(8) :: zorigin, dz_out, dy_out

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

    write(int_to_str,  "(I4.4)") outputno
    write(rank_to_str, "(I4.4)") myrank

    if (present(restart_init)) then
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
      cmd = 'mkdir -p "'//trim(subdir)//'"'
      call execute_command_line(cmd)
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
    dy_out  = dy
#else
    dy_out  = 1.0
#endif

    solnVar(1:, iGlo:, jGlo:, kGlo:) => mainVar

    ! ----------------------------
    ! Write per-rank .vti
    ! ----------------------------
    ionum = 10
    open(unit=ionum, file=vti_name, status="replace", action="write", &
         access="stream", form="unformatted")

    call write_xml(ionum, '<?xml version="1.0"?>')
    call write_xml(ionum, '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian" header_type="UInt64">')
    call write_xml(ionum, '  <ImageData WholeExtent="'//trim(extent_str(whole_i0,whole_i1,whole_j0,whole_j1,whole_k0,whole_k1))// &
                       '" Origin="'//trim(real3_str(xmin,ymin,zorigin))// &
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

    call write_xml(ionum, '      </CellData>')
    call write_xml(ionum, '    </Piece>')
    call write_xml(ionum, '  </ImageData>')
    ! binary fields follow the underscore
    call write_xml(ionum, '  <AppendedData encoding="raw">')
    write(ionum) '_'
    nbytes = 8_int64 * int(ihi-ilo+1,int64) * int(jhi-jlo+1,int64) * int(khi-klo+1,int64)
    ! VTK order: x first, then y, then z
    allocate(field(ihi-ilo+1,jhi-jlo+1,khi-klo+1))
    do ivar = NVAR_BEGIN, NVAR_END
      field = solnVar(ivar,ilo:ihi,jlo:jhi,klo:khi)
      write(ionum) nbytes
      write(ionum) field
    end do
    deallocate(field)
    write(ionum) achar(10)
    call write_xml(ionum, '  </AppendedData>')
    call write_xml(ionum, '</VTKFile>')

    close(ionum)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! ----------------------------
    ! Write .pvti on rank 0 (deterministic extents)
    ! ----------------------------
    if (myrank == 0) then
      ionum = 20
      open(unit=ionum, file=pvti_name, status="replace", action="write", form="formatted")

      write(ionum,'(A)') '<?xml version="1.0"?>'
      write(ionum,'(A)') '<VTKFile type="PImageData" version="0.1" byte_order="LittleEndian">'
      write(ionum,'(A)') '  <PImageData WholeExtent="'//trim(extent_str(whole_i0,whole_i1,whole_j0,whole_j1,whole_k0,whole_k1))// &
                         '" Origin="'//trim(real3_str(xmin,ymin,zorigin))// &
                         '" Spacing="'//trim(real3_str(dx,dy_out,dz_out))//'">'

      !write(ionum,'(A)') '      <FieldData>'
      !write(ionum,'(A)') '        <DataArray type="Int32" Name="CYCLE" NumberOfTuples="1" format="ascii">'
      !write(ionum,'(I0)') step
      !write(ionum,'(A)') '        </DataArray>'
      !write(ionum,'(A)') '        <DataArray type="Float64" Name="TIME" NumberOfTuples="1" format="ascii">'
      !write(ionum,'(ES25.18)') dble(t)
      !write(ionum,'(A)') '        </DataArray>'
      !write(ionum,'(A)') '      </FieldData>'

      write(ionum,'(A)') '    <PCellData Scalars="dens">'
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
      write(ionum,'(A)') '    </PCellData>'

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

        write(rank_to_str, "(I4.4)") r

        if (present(restart_init)) then
          write(ionum,'(A)') '    <Piece Extent="'//trim(extent_str(ilo_r-1,ihi_r, jlo_r-1,jhi_r, klo_r-1,khi_r))// &
                             '" Source="vti_reinit_'//trim(adjustl(int_to_str))//'/'//trim(adjustl(basenm))// &
                             '_reinit_'//trim(adjustl(rank_to_str))//'_'//trim(adjustl(int_to_str))//'.vti"/>'
        else
          write(ionum,'(A)') '    <Piece Extent="'//trim(extent_str(ilo_r-1,ihi_r, jlo_r-1,jhi_r, klo_r-1,khi_r))// &
                             '" Source="vti_'//trim(adjustl(int_to_str))//'/'//trim(adjustl(basenm))// &
                             '_'//trim(adjustl(rank_to_str))//'_'//trim(adjustl(int_to_str))//'.vti"/>'
        endif
      end do

      write(ionum,'(A)') '  </PImageData>'
      write(ionum,'(A)') '</VTKFile>'

      close(ionum)
    end if

    nullify(solnVar)

  end subroutine write_output


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
    write(u,'(A)') '      <PDataArray type="Float64" Name="'//trim(name)//'" NumberOfComponents="1"/>'
  end subroutine write_p_scalar

  subroutine write_scalar(u, name, solnVar, ivar, ilo, ihi, jlo, jhi, klo, khi)
    integer, intent(in) :: u, ivar, ilo, ihi, jlo, jhi, klo, khi
    character(len=*), intent(in) :: name
    real, pointer :: solnVar(:,:,:,:)
    integer(int64) :: nbytes, offset
    character(len=32) :: offset_str

    nbytes = 8_int64 * int(ihi-ilo+1,int64) * int(jhi-jlo+1,int64) * int(khi-klo+1,int64)
    offset = int(ivar-NVAR_BEGIN,int64) * (8_int64 + nbytes)
    write(offset_str,'(I0)') offset
    call write_xml(u, '        <DataArray type="Float64" Name="'//trim(name)// &
                   '" NumberOfComponents="1" format="appended" offset="'//trim(offset_str)//'"/>')
  end subroutine write_scalar

  subroutine write_xml(u, line)
    integer, intent(in) :: u
    character(len=*), intent(in) :: line
    write(u) trim(line)//achar(10)
  end subroutine write_xml

end module io_module
