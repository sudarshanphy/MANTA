module sim_data
#include "param.h"
  use mpi
  implicit none
  real(8), save :: xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz
  real(8), save :: t0, tf, dt, out_dt
  real(8), save :: gamma, cfl, grav
  real(8), save :: smalld, smallp, smalle, small
  integer, save :: nx, ny, nz, Gpts, xTpts, yTpts, zTpts, restart_no, &
                   restart_step, nend, xlbc_int, ylbc_int, xrbc_int, &
                   yrbc_int, zlbc_int, zrbc_int 
  character(len=256), save :: xlbctype, ylbctype, xrbctype, yrbctype, &
                              zlbctype, zrbctype
  character(len=256), save :: outdir, &
                              basenm, recon_method, flux_solver
  logical, save :: restart, usegrav
  real(8), parameter :: PI = 4.0 * atan(1.0)
  real(8), parameter :: smallf = 1.0e-30

  ! array which has the solution
  real(8), save, dimension(:,:,:,:), target, allocatable :: mainVar

  ! only used by MHD runs 
  real(8), save :: ch !speed for divergence cleaning


  ! MPI stuff
  integer, save :: myrank, nprocs, comm, ierr, errcode
  integer, save :: xrank, yrank, zrank
  integer, dimension(MPI_STATUS_SIZE) :: status1 

  integer, save :: ilo, ihi, jlo, jhi, klo, khi, iGlo, iGhi, jGlo, jGhi, kGhi, kGlo
  integer, save :: xblk, yblk, zblk, lnx, lny, lnz, lxTpts, lyTpts, lzTpts
  ! local block xmin and xmax
  real(8), save :: lxmin, lxmax, lymin, lymax, lzmin, lzmax
  ! blocks to the left, right, top, bottom, front, back
  integer, save :: l_blk, r_blk, t_blk, b_blk, front_blk, back_blk
  logical, save :: at_xlboundary, at_ylboundary, at_xrboundary, at_yrboundary, &
                   at_zlboundary, at_zrboundary

  integer, save :: nghost(MAXDIM)
  integer, save, dimension(MAXDIM) :: dtype_mpi
end module sim_data
