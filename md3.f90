        
        subroutine md_langevin(pos,pot,fx,fy,fz,velo)
        use global
        implicit none

        include 'md2.inc1'

        double precision      :: rand   ! random number 
!        double precision      :: gamma  ! friction coefficient, defined in common
        double precision      :: Ti    ! instance tempr
        double precision      :: gfric , noise  !

write(*,*) "In md"
        
        ! store the pos and v at t = 0
        pos0 = pos
        velo0  = velo 
        diffco = 0.0d0
        t_refer = 0.5d0

        call md_files('open')


            ! 0. initialize, done in initpos() and initvelo()
            !  calculate and print energy for initial state
                Ek = 0.0d0
                do i = 1 , totnum
                     Ek = Ek + (velo(1,i)**2 + velo(2,i)**2 + velo(3,i)**2)
                end do
                Ek = 0.5d0*Ek
                Ep = pot
                Etot = Ek + Ep

                call md_files('we')

        ctime = dt
        hdt = dt*0.5d0
        flag1 = 1   ! flag for store pos0
        !gamma = 1.0d0  ! I don't know why choose as 1, may be we can have a test
        ! gamma will be read in from input.in
        ! it is suggested gamma = 0.5/RelaxTime
        
        ! where is gamma?
        gfric = 1.0d0 - gamma*hdt
        noise = sqrt(6.0d0*gamma*tempr/dt)  ! max of random force 
        ! seems to be unphysical becasue it ~ 35, quite big force
!        noise = sqrt(6.0d0*gamma*tempr*kb)

        do while (ctime < tott ) 
            ! sample msd vcf
            if (ctime .ge. t_refer) then   !  set t @ 0.8 as the reference point
                call sample_msd_vcf(pos,velo,pos0, velo0, &
                                 diffco,vcf, msd, flag1, ctime)
            end if
        ! 1. compute a, same as f
            ! don't have to do anything!
        ! 2. compute v(t+dt/2)   2014-03-29: add friction, gfric
              do i = 1 , totnum
                  velo(1,i) = velo(1,i)*gfric + hdt*fx(i)
                  velo(2,i) = velo(2,i)*gfric + hdt*fy(i)
                  velo(3,i) = velo(3,i)*gfric + hdt*fz(i)
              end do 
            ! 3. update pos 
              do i = 1 , totnum
                   pos(:,i) = pos(:,i) + velo(:,i)*dt
!                   pos(2,i) = pos(2,i) + velo(2,i)*dt
!                   pos(3,i) = pos(3,i) + velo(3,i)*dt
                     ! chenk PBC 
                  ! if (PBC .EQ. 1) then
                    call pbc_(pos(1,i), lx, 1)  ! in common.f90
                    call pbc_(pos(2,i), ly, 1)  ! in common.f90
                    call pbc_(pos(3,i), lz, 1)  ! in common.f90
                  ! end if 
              end do 
                
    

            ! 4. update f and pot  using pos
                call  LJPE(pos,pot,fx,fy,fz)   ! in LJPE.f90

                ! add random force
                do i = 1 , totnum
                    call random_number(rand)
                    fx(i) = fx(i) + 2.0d0*noise*(rand-0.5d0)
                    call random_number(rand)
                    fy(i) = fy(i) + 2.0d0*noise*(rand-0.5d0)
                    call random_number(rand)
                    fz(i) = fz(i) + 2.0d0*noise*(rand-0.5d0)
                end do

            ! 5. update v
                Ek = 0.0d0
            do i = 1 , totnum
                velo(1,i) = velo(1,i)*gfric + hdt*fx(i)
                velo(2,i) = velo(2,i)*gfric + hdt*fy(i)
                velo(3,i) = velo(3,i)*gfric + hdt*fz(i)
                Ek = Ek + (velo(1,i)**2 + velo(2,i)**2 + velo(3,i)**2)
            end do 
                Ek = Ek*0.5d0
            ! 6. calculate and print energy
                Ep = pot
                Etot = Ek + Ep
                !calc temperature
!                Ti = Ek/totnum/1.5d0

                call md_files('we')

            ! 7. go to next step
            ctime = ctime + dt


        end do ! while (ctime < tott ) 



        ! print diffusion coefficence on the screen
        diffco = diffco*dt/3.0d0
        write(*,*) "Diffusion Coeff from VCF" , diffco/diffco2reduceddiffco 
        write(*,*) "Diffusion Coeff reduced unit" , diffco

        call md_files('close')
write(*,*) "Friction: ", gfric
write(*,*) "Max noise force", noise/force2reducedforce



        return
        end subroutine   ! md_langevin


