subroutine read_par()
        ! this subroutine reads from the parameter file and 
        ! passes the values as output
        use sim_data 
        implicit none

        ! local varaibles to read the parameter file.
        character(len=100) :: buffer, label
        integer :: pos
        integer, parameter :: no = 15 ! file unit number
        integer :: ios = 0
        integer :: line = 0 ! tracks number of lines in the paramter file

        ! set some default values, incase the user doesn't provide
        xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0; zmin = 0.0; zmax = 1.0
        nx = 10; ny = 10; nz=10; restart_no = 0; restart_step = 0
        t0 = 0.0; tf = 0.5; dt = 1.0e-4; out_dt = 1.0e-3
        gamma = 1.4; cfl = 0.4; grav = 0.0
        xlbctype = "flow"; ylbctype = "flow"; zlbctype = "flow" 
        xrbctype = "flow"; yrbctype = "flow"; zrbctype = "flow"
        nend = 10000; 
        basenm = "def"; restart = .false.; usegrav = .false.;
        recon_method = "weno3"; flux_solver = "hlle";
        xblk = 1; yblk = 1; zblk = 1; outdir = "./"
 
        open(no, file='./input.par')

        ! ios is negative if an end of record condition is encountered or if
        ! an endfile condition was detected.  It is positive if an error was
        ! detected.  ios is zero otherwise.

        do while (ios == 0)
           read(no, '(A)', iostat=ios) buffer ! read the line
           if (ios == 0) then
              line = line + 1 ! count the line

              ! skip empty lines and comments
              buffer = adjustl(buffer)
              if (len_trim(buffer) == 0) cycle
              if (buffer(1:1) == '#' .or. buffer(1:1) == '!') cycle
              
              ! split the label and value at =
              pos = scan(buffer, '=')
              if (pos <= 1) cycle
              label = trim(adjustl(buffer(:pos-1))) ! parameter name
              buffer = buffer(pos+1:) ! parameter value

              select case (label) ! based on the label, now assign values to predefined varaibles
              case ('xmin')
                 read(buffer, *, iostat=ios) xmin
              case ('xmax')
                 read(buffer, *, iostat=ios) xmax
              case ('ymin')
                 read(buffer, *, iostat=ios) ymin
              case ('ymax')
                 read(buffer, *, iostat=ios) ymax
              case ('zmin')
                 read(buffer, *, iostat=ios) zmin
              case ('zmax')
                 read(buffer, *, iostat=ios) zmax
              case ('xblk')
                 read(buffer, *, iostat=ios) xblk
              case ('yblk')
                 read(buffer, *, iostat=ios) yblk
              case ('zblk')
                 read(buffer, *, iostat=ios) zblk
              case ('nx')
                 read(buffer, *, iostat=ios) nx
              case ('ny')
                 read(buffer, *, iostat=ios) ny
              case ('nz')
                 read(buffer, *, iostat=ios) nz
              case ('t0')
                 read(buffer, *, iostat=ios) t0
              case ('tf')
                 read(buffer, *, iostat=ios) tf
              case ('dt')
                 read(buffer, *, iostat=ios) dt
              case ('out_dt')
                 read(buffer, *, iostat=ios) out_dt
              case ('gamma')
                 read(buffer, *, iostat=ios) gamma
              case ('cfl')
                 read(buffer, *, iostat=ios) cfl
              case ('grav')
                 read(buffer, *, iostat=ios) grav
              case ('usegrav')
                 read(buffer, *, iostat=ios) usegrav
              case ('xlbctype')
                 read(buffer, *, iostat=ios) xlbctype
              case ('xrbctype')
                 read(buffer, *, iostat=ios) xrbctype
              case ('ylbctype')
                 read(buffer, *, iostat=ios) ylbctype
              case ('yrbctype')
                 read(buffer, *, iostat=ios) yrbctype
              case ('zlbctype')
                 read(buffer, *, iostat=ios) zlbctype
              case ('zrbctype')
                 read(buffer, *, iostat=ios) zrbctype
              case ('nend')
                 read(buffer, *, iostat=ios) nend
              case ('basenm')
                 read(buffer, *, iostat=ios) basenm
              case ('recon_method')
                 read(buffer, *, iostat=ios) recon_method
              case ('flux_solver')
                 read(buffer, *, iostat=ios) flux_solver
              case ('restart')
                 read(buffer, *, iostat=ios) restart
              case ('restart_no')
                 read(buffer, *, iostat=ios) restart_no
              case ('restart_step')
                 read(buffer, *, iostat=ios) restart_step
              case ('output_dir')
                 read(buffer, *, iostat=ios) outdir
              case default
                ! for comments on the paramter file                
              end select
           end if
        end do
        close(no) !close the file
end subroutine
