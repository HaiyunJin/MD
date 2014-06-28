      ! LJPE 
      ! By Haiyun Jin
      ! Date:  2014-01
      ! Note: potential energy and force can be seperate into two
      ! different subroutine.  2014-03-12
      !
      !    Two subroutines:
      !  pot_s: calculation en related to an single partical, without shift
      !  LJPE : calculation total pot and force, with shiftPshiftF   
      !
      !    One function:
      !  LJPEoPBC:  LJ potential without PBC or shift
      !


     ! ! potential energy without shift, single partical
      subroutine pot_s(r,pos,no,pot,vir)  
       use global
      implicit none

      ! external variables
      integer          :: no   !  number of partical that is chosen
      double precision :: r(3), pos(3,totnum), pot, vir

      ! internal variables
      double precision :: rxij, ryij, rzij
      double precision :: dist, dist2, distr6, distr12
      integer          :: i         ! index

       call switch(pos, 1, no) ! switch no to the 1st site, in LJPE.f90

        ! reset the vir and potential
        pot = 0.0d0
        vir = 0.0d0

            do i = 2, totnum  ! ignore itself
                rxij = pos(1,i) - r(1)  
                ryij = pos(2,i) - r(2)  
                rzij = pos(3,i) - r(3) 
                ! Periodic Boundary Condition
                call pbc_(rxij, lx, 2)   ! in common.f90
                call pbc_(ryij, ly, 2)   ! in common.f90
                call pbc_(rzij, lz, 2)   ! in common.f90
                ! calc distance 
                dist2 = rxij**2 + ryij**2 + rzij**2
                !dist  = sqrt(dist2)
                distr6=1.0d0/dist2**3     ! for calc efficiency
                distr12=distr6**2
                ! calc potential energy and vir
                pot = pot +  4.0d0*(distr12 - distr6) ! no shift
                vir = vir + 16.0d0*(distr12 - distr6*0.5d0 )  ! 48/3
            end do 

       call switch(pos, 1, no) ! switch back, in LJPE.f90

        return
      end subroutine ! pot_s




        subroutine switch(pos, a, b )
        use global
        implicit none
        double precision  :: pos(3, totnum), temp(3)
        integer           :: a, b
               temp = pos(:,a)
           pos(:,a) = pos(:,b)
           pos(:,b) = temp
        return
        end subroutine  ! switch






      subroutine LJPE(pos,pot,fx,fy,fz)
       use global
      implicit none

      ! external variables
      double precision ::  pos(3,totnum), pot
      double precision ::  fx(totnum),fy(totnum),fz(totnum) 

      ! internal variables
      !double precision ::  rc, rcr6, urc, dudrc   ! cutoff of LJ curve
      double precision ::  dist, distr6, distr12
      double precision ::  rxij, ryij, rzij
      double precision ::  force, temp1, temp2
      integer   :: i, j      ! index


        ! reset the force and potential
         do i = 1, totnum
            fx(i) = 0.0d0
            fy(i) = 0.0d0
            fz(i) = 0.0d0
         end do
         pot = 0.0d0

     ! get total pot in reduced unit w/ PBC
        do i = 1, totnum
            do j = i+1, totnum
                rxij = pos(1,i)-pos(1,j)
                ryij = pos(2,i)-pos(2,j)
                rzij = pos(3,i)-pos(3,j)
                ! Periodic Boundary Condition
                call pbc_(rxij, lx, 2)   ! in common.f90
                call pbc_(ryij, ly, 2)   ! in common.f90
                call pbc_(rzij, lz, 2)   ! in common.f90

                ! calc distance 
                dist = sqrt(rxij**2 + ryij**2 + rzij**2)
                distr6=1.0d0/dist**6     ! for calc efficiency
                distr12 = distr6**2

                ! calc potential energy
                pot  = pot + 4.0d0*(distr12- distr6) - urc  &
                           - (dist - rc)*dudrc

                ! calculate force:: in reduced unit, same as accelara
                temp1 = 48.0d0/dist**2*(distr12 - distr6*0.5d0) + dudrc/dist
                force = rxij*temp1
                  fx(i) = fx(i) + force
                  fx(j) = fx(j) - force
                force = ryij*temp1
                  fy(i) = fy(i) + force
                  fy(j) = fy(j) - force
                force = rzij*temp1
                  fz(i) = fz(i) + force
                  fz(j) = fz(j) - force
            end do
        end do

        return
      end subroutine  ! LJPE()





      ! Lennard-Jones Potentail without PBC
      ! By Haiyun Jin

      function LJPEoPBC(pos)
        use global
      implicit none
      ! external variables
      double precision :: pos(3,*)
      double precision :: LJPEoPBC
      ! internal variables
      double precision :: dist , distr6
      integer :: i, j

     ! get total pot in reduced unit
        LJPEoPBC= 0.0d0
        do i = 1, totnum
            do j = i+1, totnum
                dist = sqrt((pos(1,i)-pos(1,j))**2 &
                           +(pos(2,i)-pos(2,j))**2 &
                           +(pos(3,i)-pos(3,j))**2 )
                distr6=1.0d0/dist**6
                LJPEoPBC = LJPEoPBC + 4.0d0*distr6*(distr6 - 1.0d0)
            end do
        end do

      end function


