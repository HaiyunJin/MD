        ! initialize lattice
        ! By Haiyun Jin  2014-1-22
        ! 
        ! Have 5 subroutine
        ! 1. initpos
        ! 2. randseed
        ! 3. initvelo
        ! 4. read_pos 
        ! 5. read_velo
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine initpos(pos)
        ! creat a supercell in fractional coordinates
          use global
        implicit none 

        ! external variables
        double precision    :: pos(3,*)  ! position and velocity

        ! internal variables
        integer ::  x, y, z, n, i     ! loop index

        ! Four basic atoms
        pos(1,1) = 0.0d0 
        pos(2,1) = 0.0d0 
        pos(3,1) = 0.0d0 

        pos(1,2) = 0.5d0 
        pos(2,2) = 0.5d0 
        pos(3,2) = 0.0d0 
        
        pos(1,3) = 0.5d0 
        pos(2,3) = 0.0d0 
        pos(3,3) = 0.5d0 
            
        pos(1,4) = 0.0d0 
        pos(2,4) = 0.5d0 
        pos(3,4) = 0.5d0 

        ! creating the coordinates
        n =  1     ! totalnumber
        do x = 1, nx
            do y = 1, ny
                do z = 1, nz
                    do i = 1, 4   ! 4 atoms. number of atom in one unit
                        pos(1,n) = pos(1,i) + dfloat(x-1)
                        pos(2,n) = pos(2,i) + dfloat(y-1)
                        pos(3,n) = pos(3,i) + dfloat(z-1)
                        n = n + 1 ! next array
                    end do
                end do
            end do
        end do

        do i = 1, totnum
           pos(1,i) = pos(1,i)/dfloat(nx)
           pos(2,i) = pos(2,i)/dfloat(ny)
           pos(3,i) = pos(3,i)/dfloat(nz)
        enddo 

        end subroutine   ! initpos()




      subroutine initvelo(velo)
        use global
        implicit none

       !  external variables
       double precision  :: velo(3,totnum)
        ! internal variables
       double precision  :: v0
        double precision :: qq1, qq2, s2
        integer          :: i  ! loop index

        ! evenly distribute velo
        v0 = sqrt(3.0d0*tempr)
        call randseed()   ! call personalized rand seed subroutine
        !call random_seed()   ! call system rand seed subroutine

        do i = 1, totnum
           do while ( .true. )
              call random_number(qq1)
              call random_number(qq2)
                qq1 = 2.0d0*qq1-1.0d0
                qq2 = 2.0d0*qq2-1.0d0

              s2 = qq1**2 + qq2**2
              if ( s2 < 1.0d0) then
                velo(1,i) =v0*2.0d0*sqrt(1.0d0-s2)*qq1
                velo(2,i) =v0*2.0d0*sqrt(1.0d0-s2)*qq2
                velo(3,i) =v0*(1.0d0 - 2.0d0*s2)
                exit
              end if
            end do
        end do

       end subroutine   ! initvelo()


         ! added 2014-03-14  Pi Day!
        subroutine read_pos(pos)
        use global
        implicit none
        ! external variables
        double precision     ::  pos(3,totnum)
        ! internal variables
        integer    :: i
        open (unit=3,file="POSI.txt",action="read")
        read(3,*)
        read(3,*)
        do i = 1, totnum
            read(3,*)  pos(:,i) !, pos(2,i),pos(3,i),
        end do
        close(3) ! read POSI.txt
        return
        end subroutine  ! read_pos



         ! added 2014-03-14  Pi Day!
        subroutine read_velo(velo)
        use global
        implicit none
        ! external variables
        double precision     ::  velo(3,totnum)
        ! internal variables
        integer    :: i
        open (unit=3,file="VELOI.txt",action="read")
        read(3,*)
        read(3,*)
        read(3,*)
        do i = 1, totnum
            read(3,*)  velo(:,i)  !, velo(2,i),pos(3,i),
        end do
        close(3) ! read VELOI.txt
        return
        end subroutine ! read_velo



