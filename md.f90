        !  Molecular Dynamics
        !  Haiyun
        !  I forgot the date I write it. But today is 2014-02-27




        subroutine  md(pos,pot,fx,fy,fz,velo,md_method)
        use global
        implicit none
        ! external variables
        double precision      :: pos(3,totnum), pot
        double precision      :: fx(totnum), fy(totnum), fz(totnum) 
        double precision      :: velo(3,totnum) !, vy(totnum), vz(totnum) 
        integer    ::  md_method
        
        select case (md_method)
            case (0)
                call md_nve(   pos,pot,fx,fy,fz,velo)    ! 0 nve
            case (1) 
                call md_andersen( pos,pot,fx,fy,fz,velo) ! 1 anderseb
            case (2)
                call md_berendsen(pos,pot,fx,fy,fz,velo) ! 2 berendsen
                write(*,*) "tau  ", tau
            case (3)
                call md_nh(pos,pot,fx,fy,fz,velo)        ! 3 nose-hoover
            case (4)
                call md_langevin( pos,pot,fx,fy,fz,velo) ! 4 langvevin
        end select

        return
        end subroutine   ! md



        subroutine md_nve(pos,pot,fx,fy,fz,velo)
            use global
            implicit none
              include 'md2.inc1'
              include 'md2.inc2'
              include 'md2.inc3'
            return
        end subroutine   ! md_nve()






        subroutine  sample_msd_vcf(pos,velo,pos0, velo0,diffco, vcf, msd, flag1)
        use global
        implicit none  
        ! external variables
        double precision   ::  pos(3,totnum), velo(3,totnum)
        double precision   ::  pos0(3,totnum), velo0(3,totnum)
        double precision   ::  diffco, msd, vcf
        integer            :: flag1
        ! internal variables
                   if (flag1 .eq. 1) then 
                      ! store the pos and v at t = tott/2
                      pos0 = pos
                      velo0  = velo 
                      flag1 = 0
                   end if 
                     ! Calc. and write Mean squared displacement
                    call msd_(pos, pos0, msd)    ! in msd.f90
                    call vcf_(velo, velo0, vcf)  ! in msd.f90
                     ! velocity correlation function  vcf
                    diffco = diffco + vcf
        return
        end subroutine   ! sample_msd_vcf



        subroutine md_files(act)
        use global
        implicit none
        ! external variable
        character(*)   :: act

        if (act .eq. "open") then
         open (unit=31,file="VCF.txt",action="write",status="replace")
           write(31,*) "Velocity Correlation Function  VCF "
           write(31,*) " A/ps **2"
           write(31,*) ""
           write(31,"(2A15)") "Time/ps", "<vt*v0>"
         open (unit=30,file="MSD.txt",action="write",status="replace")
           write(30,*) "Mean Squared Displacement m.s.d"
           write(30,*) " Angstrom"
           write(30,*) ""
           write(30,"(3A15)") "Time/ps", "M.S.D/A", "Diffusion Coeff"
         open (unit=2,file="ENERGY.txt",action="write",status="replace")
            write(2,"('    Total time: ',F9.5, ' ps')") tott/time2reducedtime
            write(2,"(5A15)") "Time/ps","Ek/eV","Ep/eV","Etot/eV","T/K"
            ! print Ek, Ep, and Etot

        else if (act .eq. 'we') then
           write(2,"(5F15.7)") ctime/time2reducedtime , &
                               Ek   /engtotnum_conv , &
                               Ep   /engtotnum_conv , &
                               Etot /engtotnum_conv , &
                               Ek/1.5d0/temptotnum_conv 
        else if (act .eq. 'close') then
           close(31)  ! VCF.txt
           close(30)  ! MSD.txt
           close(2)   ! ENERGY.txt
        end if 

            return
        end subroutine  ! md_files

