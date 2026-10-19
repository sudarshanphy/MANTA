module mpi_func
#include "param.h"
  implicit none
  contains
  
    subroutine init_procs()
      use sim_data, only: comm, myrank, nprocs, &
                          ierr
      use mpi 
      implicit none

      comm = MPI_COMM_WORLD ! used frequently
      ! initialize mpi
      call MPI_INIT(ierr)
      
      ! get number of processors
      call MPI_COMM_SIZE(comm, nprocs, ierr)
      
      ! get processor rank
      call MPI_COMM_RANK(comm, myrank, ierr)

    end subroutine
    
    subroutine check_procs()
      use sim_data, only: comm, myrank, nprocs, &
                          ierr, nx, ny, nz, xblk, yblk, zblk
      use mpi 
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

      ! check counts before division
      if (min(nx,ny,nz,xblk,yblk,zblk) < 1) then
        if (myrank == MASTER_PROC) print *, "Grid and block counts must be positive"
        call MPI_Abort(comm, 1, ierr)
      end if

      if (mod(nx,xblk) /= 0 .or. mod(ny,yblk) /= 0 .or. mod(nz,zblk) /= 0) then
        if (myrank == MASTER_PROC) print *, "Grid sizes must be divisible by block counts"
        call MPI_Abort(comm, 1, ierr)
      end if

      if (nprocs /= xblk * yblk * zblk) then
        if (myrank == MASTER_PROC) then
          print *, "nprocs, xblk, yblk, zblk = ", nprocs, xblk, yblk, zblk
          print *, "number of processors should be = xblk * yblk * zblk"
        end if
        call MPI_Abort(comm, 1, ierr)
      end if

    end subroutine

    subroutine finalize_procs()
      use sim_data, only: ierr, dtype_mpi
      use mpi
      implicit none
      integer :: dir

      ! release exchange types
      do dir = 1, NDIM
        if (dtype_mpi(dir) /= MPI_DATATYPE_NULL) then
          call MPI_TYPE_FREE(dtype_mpi(dir), ierr)
        end if
      end do

      ! finalize mpi
      call MPI_FINALIZE(ierr)
    end subroutine

    subroutine create_mpiDatatypes()
      use sim_data, only: ierr, dtype_mpi, nghost, &
                          lxTpts, lyTpts, lzTpts
      use mpi
      implicit none
      integer :: exch1, exch2, count
      integer(kind=MPI_ADDRESS_KIND) :: real_lb, real_extent, plane_stride

      dtype_mpi = MPI_DATATYPE_NULL

#if NDIM == 1
      ! left and right
      count = NVAR_NUMBER * nghost(IAXIS)
      call MPI_TYPE_CONTIGUOUS(count, MPI_DOUBLE_PRECISION, &
                               dtype_mpi(IAXIS), ierr)
      call MPI_TYPE_COMMIT(dtype_mpi(IAXIS), ierr)
#else
      ! x strip across y rows
      call MPI_TYPE_VECTOR(lyTpts, NVAR_NUMBER*nghost(IAXIS), &
                           NVAR_NUMBER*lxTpts, MPI_DOUBLE_PRECISION, &
                           exch1, ierr)

      ! y slab across x rows
      count = NVAR_NUMBER * lxTpts * nghost(JAXIS)
      call MPI_TYPE_CONTIGUOUS(count, MPI_DOUBLE_PRECISION, exch2, ierr)

#if NDIM == 2
      dtype_mpi(IAXIS) = exch1
      dtype_mpi(JAXIS) = exch2
      call MPI_TYPE_COMMIT(dtype_mpi(IAXIS), ierr)
      call MPI_TYPE_COMMIT(dtype_mpi(JAXIS), ierr)
#else
      ! byte stride between xy planes
      call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, real_lb, real_extent, ierr)
      plane_stride = int(NVAR_NUMBER, MPI_ADDRESS_KIND) * &
                     lxTpts * lyTpts * real_extent

      ! repeat x strips through z
      call MPI_TYPE_CREATE_HVECTOR(lzTpts, 1, plane_stride, &
                                   exch1, dtype_mpi(IAXIS), ierr)
      call MPI_TYPE_COMMIT(dtype_mpi(IAXIS), ierr)

      ! repeat y slabs through z
      call MPI_TYPE_CREATE_HVECTOR(lzTpts, 1, plane_stride, &
                                   exch2, dtype_mpi(JAXIS), ierr)
      call MPI_TYPE_COMMIT(dtype_mpi(JAXIS), ierr)

      ! z slab across xy planes
      count = NVAR_NUMBER * lxTpts * lyTpts * nghost(KAXIS)
      call MPI_TYPE_CONTIGUOUS(count, MPI_DOUBLE_PRECISION, &
                               dtype_mpi(KAXIS), ierr)
      call MPI_TYPE_COMMIT(dtype_mpi(KAXIS), ierr)

      ! release intermediate types
      call MPI_TYPE_FREE(exch1, ierr)
      call MPI_TYPE_FREE(exch2, ierr)
#endif
#endif
    end subroutine create_mpiDatatypes
end module mpi_func
