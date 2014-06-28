        !  andersen and berendsen thermostats
        !  By Haiyun 
        ! 2014-03-29
        !  include Four subroutines
        ! 1. md_andersen
        ! 2. md_berenden
        ! 3. andersen_
        ! 4. berenden_




        subroutine md_andersen(pos,pot,fx,fy,fz,velo)
        use global
        implicit none

        ! andersen variable
        double precision    :: uplim, aaa, bbb, maxvelo
        common /maxboltz/ uplim, aaa, bbb, maxvelo

        include 'md2.inc1'
        
          ! constant for maxboltz  
          uplim = sqrt(0.5d0*beta/pi)
          aaa = sqrt(0.5d0*beta/pi)
          bbb = -0.5d0*beta
          maxvelo = 10.0d0*velo2reducedvelo    ! average velo, 4

        call randseed()

        include 'md2.inc2'
        call andersen_(velo) 
        include 'md2.inc3'
write(*,*) "Update freq:",  frequenc


        return
        end subroutine   ! md_andersen



        subroutine md_berendsen(pos,pot,fx,fy,fz,velo)
        use global
        implicit none
        ! berendsen variables

        include 'md2.inc1'

        !tau = 1.0d0
        include 'md2.inc2'
        call berendsen_(velo,1)  ! 1, berendsen, 2, simple scale
        include 'md2.inc3'

write(*,*) " Coupling time:", tau
        return
        end subroutine !  md_berendsen



        subroutine andersen_(velo)
        use global
        implicit none
        !external variabeles
        double precision  :: velo(3,*)
        ! internal variables
        double precision  :: rand
        integer           ::  no   ! number of atom chosen
        double precision  ::  velon(3), maxboltz, velon2
        integer           :: i  ! index
        double precision  :: uplim, aaa, bbb, maxvelo
        common /maxboltz/ uplim, aaa, bbb, maxvelo
        

        do no = 1 , totnum
            call random_number(rand)
            if (frequenc .gt. rand) then ! scale or not
                ! 1. randomly pick a velo from 
                !    maxwell-boltzmann distribution, use reject scheme
                do i = 1 , 3
                    do while( .true. )
                        call random_number(rand)
                        rand = 2.0d0*rand - 1.0d0
                        velon(i) = rand * maxvelo
                        velon2 = velon(i)**2
                        maxboltz = aaa*exp(bbb*velon2)
                        call random_number(rand)
                        if (uplim*rand .lt. maxboltz) then
                           !write(1000,*) "Accept"
                            exit
                        end if
                        !write(1000,*) "Reject"
                    end do
                end do
                ! 2. Give new velo to selected atom
                velo(:,no) = velon(:)
            end if
        end do

        return
        end subroutine ! andersen_



        subroutine  berendsen_(velo,method)
        use global
        implicit none
        ! external variabeles
        double precision  :: velo(3,*) 
        ! internal variables
        double precision  :: lambda
        double precision  :: Ti   ! instant temperature
        integer           :: method  ! 1 is Berendson, other are simple
        integer     ::  i

 
        Ti = Ek/totnum/1.5d0  !  instant temperature
        if (method .eq. 1) then
            ! Berendsen thermostat */
            lambda = sqrt(1.0d0+dt/tau*(tempr/Ti-1.0d0))
        else
            ! Simple thermostat */
            lambda = sqrt(tempr/Ti)
        end if
        ! Rescale
        do i = 1, totnum
            velo(:,i) = velo(:,i)*lambda
        end do

        return
        end subroutine  ! berendsen_

