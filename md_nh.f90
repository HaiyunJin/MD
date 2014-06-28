
        subroutine md_nh(pos,pot,fx,fy,fz,velo)

        use global 
        implicit none
    
          include 'md2.inc1'

        double precision   ::  xi1,xi2,vxi1,vxi2
        double precision   :: dt_2, dt_4, dt_8
        common /nh/ dt_2, dt_4, dt_8

write(*,*) "In md_nh"

        dt_2 = dt/2.0d0
        dt_4 = dt_2/2.0d0
        dt_8 = dt_4/2.0d0

                
        ! store the pos and v at t = 0
        pos0 = pos
        velo0  = velo
        diffco = 0.0d0
        t_refer = 0.5d0
        
        ! initialize ksi(squigle)
        xi1 = 0.0d0
        xi2 = 0.0d0
        vxi1 = 0.0d0
        vxi2 = 0.0d0



        call md_files('open')
        
        ! initial in main.f90, now calc ener.
        do i = 1, totnum
            Ek = Ek + (velo(1,i)**2 + velo(2,i)**2 + velo(3,i)**2)
        end do
            Ek = Ek*0.5d0
            Ep = pot
            Etot = Ek + Ep
  
        ! write ener
        call md_files('we')

        ctime = dt
        flag1 = 1


        do while (ctime < tott)
!        do while (abs(ctime) < abs(tott))
            ! sample msd vcf
            if (ctime .ge. t_refer) then   !  set t @ 0.8 as the reference point
                 call sample_msd_vcf(pos,velo,pos0, velo0, &
                      diffco,vcf, msd, flag1)
            end if

        call nhchain(xi1,xi2,vxi1,vxi2,velo)

        !integrate first half time step
        do i = 1, totnum
            pos(:,i) = pos(:,i) + velo(:,i)*dt_2
            call pbc_(pos(1,i), lx, 1)  ! in common.f90
            call pbc_(pos(2,i), ly, 1)  ! in common.f90
            call pbc_(pos(3,i), lz, 1)  ! in common.f90
        end do

        ! update f at half time
        call  LJPE(pos,pot,fx,fy,fz)   ! in LJPE.f90
        
        !integrate second half time step
        Ek = 0.0d0
        do i = 1 , totnum
            ! update velo
            velo(1,i) = velo(1,i) + dt*fx(i)
            velo(2,i) = velo(2,i) + dt*fy(i)
            velo(3,i) = velo(3,i) + dt*fz(i)
            ! update pos
            pos(:,i) = pos(:,i) + velo(:,i)*dt_2
            ! pbc will be call in the next step
            !call pbc_(pos(1,i), lx, 1)  ! in common.f90
            !call pbc_(pos(2,i), ly, 1)  ! in common.f90
            !call pbc_(pos(3,i), lz, 1)  ! in common.f90
            ! calc knetic energy
            Ek = Ek + (velo(1,i)**2 + velo(2,i)**2 + velo(3,i)**2)
        end do
            Ek = Ek*0.5d0
            Ep = pot

        call nhchain(xi1,xi2,vxi1,vxi2,velo)

            Etot = Ek + Ep
  
            ! print current ener
            call md_files('we')
 
          ctime = ctime + dt  ! go to next step
        end do ! while(t<tott)


        call md_files("close")

write(*,*) "NH mass", Q1, Q2


        return
        end subroutine   ! md_nh


        subroutine nhchain(xi1,xi2,vxi1,vxi2,velo)
        use global
        implicit none

        ! external variables
        double precision   :: xi1,xi2,vxi1,vxi2
        double precision   :: velo(3,totnum)
        ! internal variables
        double precision   :: G1, G2, s
        double precision   :: dt_2, dt_4, dt_8
        common /nh/ dt_2, dt_4, dt_8
        integer  :: i ! index

        G2 = Q1*vxi1*vxi1 - tempr
        vxi2 = vxi2 + G2 * dt_4  ! error
        vxi1 = vxi1 * exp(-vxi2*dt_8)
        G1 = (2.0d0*Ek - 3.0d0*totnum*tempr)/Q1
        vxi1 = vxi1 + G1*dt_4
        vxi1 = vxi1*exp(-vxi2*dt_8)
        xi1 = xi1 + vxi1*dt_2
        xi2 = xi2 + vxi2*dt_2
        s = exp(-vxi1*dt_2)
        do i = 1, totnum
            velo(:,i) = velo(:,i)*s
        end do
        Ek = Ek*s*s
        vxi1 = vxi1 * exp(-vxi2*dt_8)
        G1 = (2.0d0*Ek - 3.0d0*totnum*tempr)/Q1
        vxi1 = vxi1 + G1*dt_4
        vxi1 = vxi1*exp(-vxi2*dt_8)
        G2 = (Q1*vxi1*vxi1 - tempr)/Q2
        vxi2 = vxi2 + G2 * dt_4

        return
        end subroutine ! nhchain 


