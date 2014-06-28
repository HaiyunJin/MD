        !   Monte Carlo 
        !   By Haiyun
        !   2014-3-12
        !
        !  Contains:
        !   mc() 
        !   imove()
        !   open_file_mc
        !
        !

        subroutine mc(pos,pot,ncycle)
        use global
        implicit none

        ! external varibles
        double precision      :: pos(3,totnum) , pot
        integer               ::  ncycle

        ! internal varibles
        double precision      :: mcpot, mcpot2
        double precision      :: pressure, mcpres         !  system pressure  
        double precision      :: cv               !  heat capacity
!        double precision      :: ! Ek  !  kinetic energy, defined in global
        double precision      :: beta_T, tot_engconv, tot_cvconv
        logical               :: ichange   ! accept or not
        integer               :: i


        call open_file_mc(ncycle)    ! in mc.f90

        call pressure_(pos,pressure)  ! in pressure.f90

        mcpot = 0.0d0
        mcpot2= 0.0d0
        mcpres= 0.0d0
        Ek    = 1.5d0*dfloat(totnum)/beta
        beta_T= beta/tempr
        tot_engconv = dfloat(totnum)*eng2reducedeng
        tot_cvconv  = tot_engconv/temp2reducedtemp

        do i = 1, ncycle
            call imove(pos,pot,pressure,ichange)   ! in mc.f90

            mcpot = (    pot  - mcpot )/dfloat(i) +mcpot
            mcpot2= (  pot**2 - mcpot2)/dfloat(i) +mcpot2
            mcpres= (pressure - mcpres)/dfloat(i) +mcpres
! After derivation, I found I don't need to include Ek, it's constant.

            cv    = (mcpot2 - mcpot**2)*beta_T

            if ((i .eq. 1 ) .or. (mod(i,100) .eq. 0 )) then
            !    write(33,'(I10,F15.7)') i,   mcpot /tot_engconv
                write(35,'(I10,F15.7,F15.7)') i, pot/tot_engconv,  (pot+Ek)/tot_engconv
                write(34,'(I10,E15.7)')   i,     cv  /tot_cvconv
                write(36,'(I10,E15.7,E15.7)') i,pressure/press2reducedpress, &
                                                  mcpres/press2reducedpress 
            end if 
        end do
        close(33)
        close(34)
        close(35)
        close(36)


        return
        end subroutine  ! mc










        subroutine imove(pos,pot,pressure,ichange)
        use global
        implicit none

        ! external varibles
        double precision      :: pos(3,totnum), pot, pressure
        logical               :: ichange   ! accept or not

        ! internal varibels
        double precision      :: opos(3), npos(3)
        double precision      :: olden, newen    ! old and new energy
        double precision      :: oldvir, newvir    ! old and new vir
        double precision      :: rand
        integer               :: no  ! randomly picked number of partical
        integer               :: i

        call randseed()     ! in common.f90
!        deltr = 0.5d0       ! choose as radius of the partical

        call random_number(rand)  ! system
        no = int( rand * totnum ) + 1    ! randomly pick one partical
        opos = pos(:,no)

        ! calculate old and new energy
        call pot_s(opos, pos, no, olden, oldvir)    ! in LJPE.f90
        do i = 1,3
            call random_number(rand)   ! system
            npos(i) = opos(i) + (rand - 0.5d0)*deltr   
        end do
          call pbc_(npos(1),lx,1)    ! PBC,  boundary ! in common.f90
          call pbc_(npos(2),ly,1)    ! PBC,  boundary ! in common.f90
          call pbc_(npos(3),lz,1)    ! PBC,  boundary ! in common.f90
        call pot_s(npos, pos, no, newen, newvir)    ! in LJPE.f90

        call random_number(rand)  ! system
        ichange = (rand .lt. exp(-beta*(newen-olden)) ) 
        if (ichange) then  ! if ichange is true, use new pos
            pos(1,no) = npos(1)
            pos(2,no) = npos(2)
            pos(3,no) = npos(3)
            pot = pot - olden  + newen  ! get the new energy
            pressure = pressure + (newvir - oldvir)/totvolm
        end if      

        return
        end subroutine! imove()





        subroutine open_file_mc(ncycle)
        implicit none
        integer               :: ncycle

        open (unit=33,file="MC_POT.txt",action="write",status="replace")
        write(33,*) "Monte Carlo"
        write(33,*) ncycle
        write(33,*) " ncycle  Potential Eng"

        open (unit=35,file="MC_ENERGY.txt",  &
                            action="write",status="replace")
        write(35,*) "Monte Carlo Total energy"
        write(35,*) ncycle
        write(35,*) " ncycle      Etot"

        open (unit=34,file="MC_CV.txt",action="write",status="replace")
        write(34,*) "Monte Carlo heat capacity"
        write(34,*) ncycle
        write(34,*) " ncycle      Cv"

        open (unit=36,file="MC_PRES.txt",action="write",status="replace")
        write(36,*) "Monte Carlo System Pressure, in eV/A^3"
        write(36,*) ncycle
        write(36,*) " ncycle        Pressure          mcpres       "


        return
        end subroutine ! open_file_mc()
