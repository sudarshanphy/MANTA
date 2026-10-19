program hydro
#include "param.h"

  use sim_data
  use grid_func
  use sim_init, only: init_problem
  use io_restart, only: restart_problem
  use io_module, only: write_output
  use applyBC_module, only: applyBC_all
  use rk2_module, only: RK2_SSP
  use misc_module, only: get_dt
#ifdef MHD
  use misc_module, only: compute_maxdivB
#endif
  use mpi_func
  use guard_func               
  use mpi 

  implicit none

  real :: time, timeio
  integer ::  step, outputno
  logical :: io_output
  
  real :: timer_start, timer_stop, time_per_step, timer_step0, &
          timer_step1
 
#ifdef MHD
  real(8) :: maxdivB, maxdivB_error
#endif

  ! initilize mpi
  call init_procs()

  !read the parameter file by each processor
  call read_par()

  !check processor numbers
  call check_procs()

  !initialize the grid
  call grid_init()
  
  call create_mpiDatatypes()
  ! check if it is a restart or not
  if (.not. restart) then
    step = 0
    outputno = 0

    !initialize the problem
    call init_problem()
                                            
    call guardcell_fill()
    call applyBC_all()
    call write_output(t0, step, outputno) 

  else

    ! restart a problem from an output file
    call restart_problem(restart_no, t0)
    step = restart_step
    outputno = restart_no
    call guardcell_fill()
    call applyBC_all()
    call write_output(t0, step, outputno, .true.) 

  end if

  call CPU_TIME(timer_start)
  timer_step0 = 0
  
  time = t0
  timeio = time + out_dt
  ! evolution loop
  do while ((time < tf) .and. (step < nend))

    io_output = .false.
    
    ! compute dt
    call get_dt(dt)

    if (time + dt > tf) then
      dt = tf - time
      io_output = .true.
    end if
    
    if (time + dt > timeio) then
      dt = timeio - time
      timeio = timeio + out_dt
      io_output = .true.
    end if


#ifdef MHD
    ! glm update psi 
    call glm(dt)
#endif 
    ! do a SSP RK2 step
    call RK2_SSP(dt)

    time = time + dt
    step = step + 1

    if (step == nend) io_output=.true.

#ifdef MHD
    ! global maxima after RK2 fills guard cells
    maxdivB = compute_maxdivB(mainVar(BMFX_VAR,:,:,:), mainVar(BMFY_VAR,:,:,:), mainVar(BMFZ_VAR,:,:,:))
    maxdivB_error = compute_maxdivB(mainVar(BMFX_VAR,:,:,:), mainVar(BMFY_VAR,:,:,:), &
                                  mainVar(BMFZ_VAR,:,:,:), normalized=.true.)
#endif
    if (myrank == MASTER_PROC) then
      print *, "Step: ", step," --> time = ", time , ", dt = ", dt
#ifdef MHD
      print '(A,ES14.6,A,ES14.6)', '  max |divB| = ', maxdivB, '  max divB_error = ', maxdivB_error
#endif
    end if

    if (io_output) then

      call CPU_TIME(timer_stop)
      timer_step1 = step
      time_per_step = (timer_stop - timer_start) &
                     /(timer_step1 - timer_step0)
      
      if (myrank == MASTER_PROC) then
        print *, "# step0 = ", timer_step0, " step1 = ", timer_step1
        print *, "# time per step = ", time_per_step
      end if
      timer_step0 = step

      outputno = outputno + 1
      call write_output(time, step, outputno)
      if (myrank == MASTER_PROC) then
        print *, "$$$ output basenm = ", trim(basenm)
        print *, "$$$ Step: ", step," || time = ", time, " || output # = ", outputno," $$$" 
      end if
      call CPU_TIME(timer_start)
    end if
#ifdef MHD
    call glm(dt)
#endif
  end do


  call grid_finalize()
  call finalize_procs()
end program hydro
