!           Mean Squared Displacement  M.S.D
!           Haiyun  
!           2014-02-26
!
!

        subroutine msd_(pos, pos0, msd)
        use global
        implicit none

        ! external variables
        double precision    :: pos(3,*), pos0(3,*), msd

        ! internal variables
        double precision    :: rx, ry, rz
        integer             :: i 

        msd = 0.0d0
        do i = 1, totnum
            rx = pos(1,i)-pos0(1,i)
            ry = pos(2,i)-pos0(2,i)
            rz = pos(3,i)-pos0(3,i)
            ! Periodic Boundary Condition, 1 for position, 2 for vector
            call pbc_(rx, lx, 2)   ! in common.f90
            call pbc_(ry, ly, 2)   ! in common.f90
            call pbc_(rz, lz, 2)   ! in common.f90

            ! calc distance 
            msd = msd + rx**2 + ry**2 + rz**2
        end do 

        msd = msd/ totnum

        write(30,"(3F15.7)") ctime/time2reducedtime, &
                             msd/length2reducedlength2, &
                             msd/ctime/6.0d0/msd_time_conv

        return
        end subroutine  ! msd_()


        !!!
        !!!
        !!!   Subroutine vcf_  <v(0)*v(t)>

        subroutine  vcf_(velo, velo0, vcf)
        use global
        implicit none

        ! external variables
        double precision    :: velo(3,totnum)   !, vy(totnum), vz(totnum)
        double precision    :: velo0(3,totnum)  !,vy0(totnum),vz0(totnum)
        double precision    :: vcf

        ! internal variables
        integer             :: i   ! loop index

        vcf = 0.0d0
        do i = 1, totnum
            vcf = vcf + velo(1,i)*velo0(1,i) &
                      + velo(2,i)*velo0(2,i) &
                      + velo(3,i)*velo0(3,i) 
        end do
        vcf = vcf /totnum
        write(31,"(2F15.7)") (ctime)/time2reducedtime, &
                                 vcf/velo2reducedvelo2

        return
        end subroutine  ! vcf_()

