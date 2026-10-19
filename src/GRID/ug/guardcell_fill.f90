module guard_func
#include "param.h"
  implicit none
  contains

  subroutine guardcell_fill()
    use sim_data, only: ilo, ihi, jlo, jhi, klo, khi, &
                        iGlo, jGlo, kGlo, mainVar, myrank, &
                        xblk, yblk, zblk, comm, ierr, status1, &
                        at_xlboundary, at_xrboundary, &
                        at_ylboundary, at_yrboundary, &
                        at_zlboundary, at_zrboundary, &
                        l_blk, r_blk, b_blk, t_blk, back_blk, front_blk, &
                        nghost, dtype_mpi
    use mpi
    implicit none
    integer :: guard
    integer, parameter :: tag1=20, tag2=21, tag3=23, tag4=24, tag5=25, tag6=26

    ! each block must supply a full guard layer
    if (ihi-ilo+1 < nghost(IAXIS) .or. &
        jhi-jlo+1 < nghost(JAXIS) .or. &
        khi-klo+1 < nghost(KAXIS)) then
      print *, "Local block is smaller than guard width on rank", myrank
      call MPI_Abort(comm, 1, ierr)
    end if

    ! inactive neighbors
    l_blk = MPI_PROC_NULL; r_blk = MPI_PROC_NULL
    b_blk = MPI_PROC_NULL; t_blk = MPI_PROC_NULL
    back_blk = MPI_PROC_NULL; front_blk = MPI_PROC_NULL

    ! wrap neighbors; physical BCs are applied afterwards
    l_blk = merge(myrank+(xblk-1), myrank-1, at_xlboundary)
    r_blk = merge(myrank-(xblk-1), myrank+1, at_xrboundary)

    ! x exchange
    guard = nghost(IAXIS)
    call MPI_SENDRECV(mainVar(1,ilo,jGlo,kGlo), 1, dtype_mpi(IAXIS), l_blk, tag1, &
                      mainVar(1,ihi+1,jGlo,kGlo), 1, dtype_mpi(IAXIS), r_blk, tag1, &
                      comm, status1, ierr)
    call MPI_SENDRECV(mainVar(1,ihi-guard+1,jGlo,kGlo), 1, dtype_mpi(IAXIS), r_blk, tag2, &
                      mainVar(1,iGlo,jGlo,kGlo), 1, dtype_mpi(IAXIS), l_blk, tag2, &
                      comm, status1, ierr)

#if NDIM >= 2
    b_blk = merge(myrank+xblk*(yblk-1), myrank-xblk, at_ylboundary)
    t_blk = merge(myrank-xblk*(yblk-1), myrank+xblk, at_yrboundary)

    ! y exchange
    guard = nghost(JAXIS)
    call MPI_SENDRECV(mainVar(1,iGlo,jlo,kGlo), 1, dtype_mpi(JAXIS), b_blk, tag3, &
                      mainVar(1,iGlo,jhi+1,kGlo), 1, dtype_mpi(JAXIS), t_blk, tag3, &
                      comm, status1, ierr)
    call MPI_SENDRECV(mainVar(1,iGlo,jhi-guard+1,kGlo), 1, dtype_mpi(JAXIS), t_blk, tag4, &
                      mainVar(1,iGlo,jGlo,kGlo), 1, dtype_mpi(JAXIS), b_blk, tag4, &
                      comm, status1, ierr)
#endif

#if NDIM == 3
    back_blk = merge(myrank+xblk*yblk*(zblk-1), myrank-xblk*yblk, at_zlboundary)
    front_blk = merge(myrank-xblk*yblk*(zblk-1), myrank+xblk*yblk, at_zrboundary)

    ! z exchange
    guard = nghost(KAXIS)
    call MPI_SENDRECV(mainVar(1,iGlo,jGlo,klo), 1, dtype_mpi(KAXIS), back_blk, tag5, &
                      mainVar(1,iGlo,jGlo,khi+1), 1, dtype_mpi(KAXIS), front_blk, tag5, &
                      comm, status1, ierr)
    call MPI_SENDRECV(mainVar(1,iGlo,jGlo,khi-guard+1), 1, dtype_mpi(KAXIS), front_blk, tag6, &
                      mainVar(1,iGlo,jGlo,kGlo), 1, dtype_mpi(KAXIS), back_blk, tag6, &
                      comm, status1, ierr)
#endif
  end subroutine guardcell_fill
end module guard_func
