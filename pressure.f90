        !   calculate pressure of the system
        !   Haiyun Jin
        !   2014-03-20

        
        subroutine  pressure_(pos,pressure)
        use global 
        implicit none

      ! external variables
      double precision ::  pos(3,totnum), pressure

      ! internal variables
      double precision ::  dist, dist2, distr6, distr12
      double precision ::  rxij, ryij, rzij
        integer    :: i , j


      ! reset pressure
        pressure = 0.0d0

     ! get pressure in reduced unit w/ PBC
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
                dist2 = (rxij**2 + ryij**2 + rzij**2)
                dist = sqrt(dist2)
                distr6=1.0d0/dist2**3     ! for calc efficiency
                distr12 = distr6**2

                if (dist .lt. rc )  then
                    !pressure = pressure - (distr12 - distr6*0.5d0)
                    pressure = pressure + (distr12 - distr6*0.5d0)
                end if 

            end do
        end do
        pressure = 16.0d0*pressure/totvolm

        pressure = pressure & 
           + 16.0d0/3.0d0 *pi*rho**2 *(2.0d0/3.0d0/rc**9 - 1.0d0/rc**3)
        write(*,*)  "Pressure - ρkBT    ", (pressure- rho/beta)/press2reducedpress 
        pressure = pressure + rho/beta
        write(*,*)  "rho*kB*T      ", rho/beta/press2reducedpress 

        return
        end subroutine !  pressure 
