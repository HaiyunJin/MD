
        ! module global
        module global
        implicit none
        double precision  :: eng2reducedeng
        double precision  :: length2reducedlength
        double precision  :: force2reducedforce
        double precision  :: temp2reducedtemp
        double precision  :: time2reducedtime
        double precision  :: velo2reducedvelo
        double precision  :: velo2reducedvelo2 ! velo square
        double precision  :: length2reducedlength2  ! square
        double precision  :: press2reducedpress
        double precision  :: diffco2reduceddiffco
        double precision  :: msd_time_conv
        double precision  :: temptotnum_conv
        double precision  :: engtotnum_conv

        double precision  :: pi, e
        parameter        ( pi = 3.141592654d0   )
        parameter        (  e = 2.718281828459d0)
        double precision  ::  lx, ly, lz, totvolm
        double precision  ::  rc, urc, dudrc
        double precision  ::  tempr , rho  ! temperature and density
        double precision  ::  kb
        double precision  ::  sigma
        double precision  ::  smallepsilon
        double precision  ::  mass
        double precision  ::  celllength
        double precision  ::  beta      !  1/kT
        double precision  ::  deltr     !  delt r in MC 
        integer ::  totnum
        integer ::  PBC   ! Periodic boundary condition
        integer ::  nx, ny, nz        ! user defined cell number
        double precision  :: dt, tott  !   time step and total running time
        double precision  :: frequenc 
        double precision  :: tau     ! tau:"rise time", time to rise to tempr

        common /andersen/ frequenc
        common /berendsen/ tau
        common /group1/ totnum, nx, ny, nz
        common /group2/ lx, ly, lz, totvolm, tempr , rho
        common /cutoff/ rc, urc, dudrc
        common /unitconvert/     eng2reducedeng, &
                           length2reducedlength, & 
                             force2reducedforce, &
                               temp2reducedtemp, &
                               time2reducedtime, &
                               velo2reducedvelo, &
                              velo2reducedvelo2, &
                          length2reducedlength2, &
                             press2reducedpress, &  
                           diffco2reduceddiffco, &
                                  msd_time_conv, & 
                                temptotnum_conv, &
                                 engtotnum_conv

        common /mc_time/     dt, tott

        double precision   :: Q1, Q2  ! nose-hoover masses
        common /nh_mass/  Q1, Q2
!        double precision    :: t_refer
!        common /md_collect/  t_refer

        double precision   :: ctime, Ek, Ep, Etot
        common /md_para/  ctime, Ek, Ep, Etot


        double precision   :: gamma  ! langevin, friction coeffic
        common /langevin/  gamma

        end module





        subroutine  set_cutoff
        use global
        implicit none
        double precision   :: rcr6

       ! rc is the smallest l/2
       rc = lx/2.0d0
       if (ly/2.0d0 < rc) then 
            rc = ly/2.0d0
       end if
       if (lz/2.0d0 < rc) then 
            rc = lz/2.0d0
       end if 
     
        ! implement shift-pot shift-force
       rcr6  =   1.0d0/rc**6     
       urc   =   4.0d0*rcr6 *(rcr6 - 1.0d0)
       dudrc = -48.0d0/rc*rcr6 *(rcr6 - 0.5d0)

        end subroutine  ! set_cutoff



          subroutine randseed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms
          
            call random_seed(size = n)
            allocate(seed(n))
            !un = 100
            ! First try if the OS provides a random number generator
            open(unit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old",iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID
               ! is
               ! useful in case one launches multiple instances of the
               ! same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
          end subroutine randseed
            ! init_random_seed






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Periodic Boundary Condition
        subroutine pbc_(a, b , c)
        double precision  ::  a, b
        integer           ::  c
        double precision  ::  left, right
        
        ! two kinds of PBC, 1 is pos, the other is vector length
        if (c .eq. 1) then
            right = b 
            left  = 0.0d0
        else
            right =  b/2.0d0
            left  = -b/2.0d0
        end if
    
        ! PBC
        if (a > right) then
            a = a - b
        else if ( a < left) then
            a = a + b
        end if

        return
        end subroutine  ! pbc_()


