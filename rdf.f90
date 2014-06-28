        ! compute g(r)
        !  Haiyun 
        !  2014-02-25     
        
        subroutine rdf(pos,filename)
        use global
        implicit none

      ! external variables
        double precision              :: pos(3,*)
        character(*)                  :: filename        

      ! internal variables
        integer                       ::  nbin  ! number of bins
        double precision, allocatable ::  gr(:)
        double precision              ::  dr, dr3    ! window width
        integer                       ::  i , j , k
        double precision              ::  rxij, ryij, rzij, dist
        double precision              ::  dvolume, nbulk

        nbin = 200
        allocate( gr(nbin) ) 


        ! get the window width
        dr = rc / dfloat(nbin)

         ! initilize 
        do k = 1, nbin
            gr(k) = 0.0d0
        end do 

        ! calculate
        do i = 1, totnum 
            do j = i+1, totnum
                rxij = pos(1,i)-pos(1,j)
                ryij = pos(2,i)-pos(2,j)
                rzij = pos(3,i)-pos(3,j)
                ! Periodic Boundary Condition
                call pbc_(rxij, lx, 2)  ! in common.f90
                call pbc_(ryij, ly, 2)  ! in common.f90
                call pbc_(rzij, lz, 2)  ! in common.f90

                ! calc distance 
                dist = sqrt(rxij**2 + ryij**2 + rzij**2)
                k  = idint(dist/dr) + 1     ! the first window is 1   
                if (k .le. nbin ) then      ! discard atoms that farther away
                    gr(k) = gr(k) + 2.0d0   ! contribution from i->j and j->i
                end if

            end do
        end do

        ! normalize gr
write(*,*)  "rho in rdf" , rho
        dr3 = dr**3
        do k = 1 , nbin
            !dvolume = (k**3 - (k-1)**3)*dr**3  ! from the book
            dvolume = 4.0d0/3.0d0 * pi * (3.0d0*dfloat(k)**2+3.0d0*dfloat(k)+1.0d0) * dr3    ! what I think.
            nbulk = dvolume * rho
            gr(k) = gr(k)/nbulk/totnum
        end do


            open(unit=29,file=filename,action="write",status="replace")
                write (29,*) "Radial Distribution Funtion"
                write (29,*)  nbin
                write (29,*)  ""
                do k = 1, nbin
                    write (29,*) dfloat(k)*dr/length2reducedlength, gr(k)
                end do  
            close(29)   ! RDF1.txt

        return
        end subroutine   ! rdf()

