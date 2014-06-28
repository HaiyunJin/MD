     
        ! By Haiyun Jin 2014-1-22        

        program main

        use global
        implicit none

        double precision, allocatable :: pos(:,:)
!        double precision, allocatable :: pos0(:,:)
        double precision, allocatable :: fx(:), fy(:), fz(:) 
        double precision, allocatable :: velo(:,:) !vx(:), vy(:), vz(:) 
        double precision              :: pot     ! total potential energy
        double precision, external    :: LJPEoPBC
        double precision              :: Ecrystal
        integer                       :: ncycle   ! cycle of mc
        integer                       :: i, j, n ! loop index
        double precision              :: scalefactor
        double precision              :: rho0
        integer                       :: md_method
        integer                       :: iniposflag, inivelflag, scaleflag
        integer           :: NumOfAtomInOneUnit
        


        !      time     in ps (SI unit)
        !      energy   in eV
        !      distance in Angst
        !      tempr    in K
        !      Thus, mass has not choice.
        !      1 unit = 1.602176487**-23kg
        !      mass = 6.6**-23 g --> 0.004119396 unit 

        ! basic value to define unit  For Argon Ar
        kb = 0.000086173324d0   ! eV/K
        sigma        = 3.4d0    ! angstrom
        smallepsilon = 0.0104d0 ! eV
        mass         = 0.004119396d0    !in eV*ps2/A2  1 unit = 1.602176487**-23kg
        NumOfAtomInOneUnit = 4
        scalefactor  = 1.0d0


        celllength   = 5.26d0   ! Angstrom

call when
        ! set up conversion
        eng2reducedeng        = 1.0d0/smallepsilon     ! from eV
        length2reducedlength  = 1.0d0/sigma           ! from angstrom
        force2reducedforce    = sigma/smallepsilon   ! eV/Angs
        temp2reducedtemp      = kb/smallepsilon
        time2reducedtime      = sqrt(smallepsilon/mass)/sigma
        velo2reducedvelo      = length2reducedlength/time2reducedtime
        velo2reducedvelo2     = velo2reducedvelo**2
        length2reducedlength2 = length2reducedlength**2
        press2reducedpress    = force2reducedforce/length2reducedlength2
        diffco2reduceddiffco  = time2reducedtime*velo2reducedvelo2

        !read parameter from input.in
        open (unit=1,file="input.in",action="read")
            read(1,*) nx, ny, nz , PBC
            read(1,*) tott, dt      ! read  in LJ unit
            read(1,*) tempr , ncycle       ! read in temperature, and cycle of mc
            read(1,*) scaleflag  ! read in scalefactor
            read(1,*)  deltr      !for mc
            read(1,*) rho0  ! desired rho of the system
            read(1,*) iniposflag , inivelflag ! initilize or read structure: 0 init, others read
            read(1,*) md_method      !for md
            read(1,*) frequenc       ! 1.frequence for andersen
            read(1,*) tau            ! 2.coupling timescale of berendsen
            read(1,*) Q1, Q2         ! 3.nose hover masses
            read(1,*) gamma          ! 4.langevin friction coeffi
        close(1)  ! read input.in
          

          ! convert to reduced unit
        celllength = celllength * length2reducedlength
             tempr =      tempr * temp2reducedtemp
              tott =       tott !* time2reducedtime
              beta = 1.0d0/(kb*eng2reducedeng/temp2reducedtemp*tempr) ! kb has the unit of energy, 
            totnum = NumOfAtomInOneUnit*nx*ny*nz
                lx = nx * celllength  
                ly = ny * celllength  
                lz = nz * celllength  
           totvolm = lx*ly*lz  ! total volume
               rho = dfloat(totnum) / totvolm  ! density 
write(*,*) "beta", beta
write(*,*) "total time", tott/time2reducedtime


        allocate( pos(3,totnum))
        allocate( fx(totnum))
        allocate( fy(totnum))
        allocate( fz(totnum))
        allocate( velo(3,totnum))
            ! additional conversion
                msd_time_conv = length2reducedlength2/time2reducedtime
              temptotnum_conv = temp2reducedtemp*totnum
               engtotnum_conv = eng2reducedeng*totnum  
!      if (abs(tott) .gt. 0) then   ! if md, then allocate
!          allocate( pos0(3,totnum))
!      end if 

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Read or Inintialize system           !!!!!!!!!!!!!!!!!!!!!       
        if (iniposflag .eq. 1 ) then 
            call read_pos(pos)   ! in initialize.f90
        else
            call initpos(pos)    ! creat the lattice 
        end if
        if (inivelflag .eq. 1 ) then 
            call read_velo(velo) ! 
        else
            call initvelo(velo)  ! initialize velocity
        end if

           ! system initial density
        if (scaleflag .eq. 1 ) then 
           scalefactor = (rho/rho0)**(1.0d0/3.0d0)
write(*,*) "scalefactor", scalefactor
           lx = lx * scalefactor 
           ly = ly * scalefactor 
           lz = lz * scalefactor 
           totvolm = lx*ly*lz  ! total volume after rescaling
           rho = dfloat(totnum) / (totvolm)  ! density after rescaling
        end if 

        ! convert fractional to reduced unit
        do i = 1 , totnum
          pos(1,i) = pos(1,i) * lx
          pos(2,i) = pos(2,i) * ly
          pos(3,i) = pos(3,i) * lz
        end do 
!!!   End of read or initialize system   !!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call set_cutoff   ! set cutoff rc.  in common.f90
write(*,*) "    rc   ", rc



        if (PBC .EQ. 1 ) then !!!!!!!!!!

           ! calculate LJ potentail energy w/ PBC
           call LJPE(pos,pot,fx,fy,fz)

Ecrystal = 0.0d0  ! calc total x direction force
do i = 1, totnum
    Ecrystal = Ecrystal + fx(i)
end do 
write(*,*) "fx total = ", Ecrystal

Ecrystal = 0.0d0  ! calc system temperature
do i = 1, totnum
  Ecrystal = Ecrystal+ (velo(1,i)**2 + velo(2,i)**2 +velo(3,i)**2)
end do 
  Ecrystal = Ecrystal /totnum
  write(*,*) "Mean v:  " , sqrt(Ecrystal), " A/ps"
  write(*,*) "System Temperature: ",  Ecrystal/3.0d0/temp2reducedtemp, "K"

Ecrystal = 0.0d0  ! calc net x motion
do i = 1, totnum
    Ecrystal = Ecrystal + velo(1,i)
end do 
write(*,*) "vx total = ", Ecrystal/totnum

Ecrystal = 0.0d0  ! y 
do i = 1, totnum
    Ecrystal = Ecrystal + velo(2,i)
end do 
write(*,*) "vy total = ", Ecrystal/totnum

Ecrystal = 0.0d0  ! z
do i = 1, totnum
    Ecrystal = Ecrystal + velo(3,i)
end do 
write(*,*) "vz total = ", Ecrystal/totnum


        else              !!!!!!!!!!
            ! calculate LJ potentail energy w/o PBC 
            pot = LJPEoPBC(pos)
        end if            !!!!!!!!!!
            Ecrystal = pot/ totnum

!        write(*,"(A,F10.6,A)") "Cohesive energy of the crystal is " , Ecrystal, " eV"
        write(*,"(3I5,F15.7,A)") nx, ny,nz,Ecrystal/eng2reducedeng," eV"





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Print structures   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! print position to POS.txt
        open (unit=100,file="POS.txt",action="write",status="replace")
        write (100,*)  totnum 
        write (100,*)  "Angstrom" 
        do n = 1, totnum
            write (100, "(3F15.7)" ) pos(:,n)/length2reducedlength
        end do
        close(100)  ! POS.txt


      ! print force of first 50 atom to FORCE.txt
         open (unit=27,file="FORCE.txt",action="write",status="replace")
         write (27,*)  "Force of first 50 atom" 
         write (27,*)  "" 
         write (27,*)  "  eV/A" 
            do i = 1, 50
               write (27, "(3F15.7)" ) fx(i)/force2reducedforce, &
                                       fy(i)/force2reducedforce, &
                                       fz(i)/force2reducedforce  
            end do
        close(27)  ! FORCE.txt


      ! print velocity to VELO.txt
!         open (unit=28,file="VELO.txt",action="write",status="replace")
!         write (28,*)  "Velocity of atoms" 
!         write (28,*)  "" 
!         write (28,*)  "  A/ps (1 A/ps = 100 m/s)" 
!            do i = 1, totnum
!               write (28, "(3F15.7)" ) velo(:,i)/velo2reducedvelo !, &
!            end do
!        close(28)  ! VELO.txt
!
!!!    End of Print structures     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       g(r) before md
!        call rdf(pos,"RDF1.txt")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           MD
        if (abs(tott) .gt. 0.0 ) then
            call md(pos,pot,fx,fy,fz,velo, md_method)
!           g(r) after md
            call rdf(pos,"RDFmd.txt")
write(*,*) " "
write(*,*) "After MD"
Ecrystal = 0.0d0  ! system tem
do i = 1, totnum
 Ecrystal = Ecrystal+ (velo(1,i)**2 + velo(2,i)**2 +velo(3,i)**2)
end do 
 Ecrystal = Ecrystal /totnum
 write(*,*) "Mean v:  " , sqrt(Ecrystal), " A/ps"
 write(*,*) "System Temp aft md: ",  Ecrystal/3.0d0/temp2reducedtemp, "K"
        end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Monte Carlo
        if (ncycle .gt. 0 ) then
            call mc(pos,pot,ncycle)  ! need to have ncycle as input
!           g(r) after mc
            call rdf(pos,"RDFmc.txt")
write(*,*) "After MC"
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! print final position to POSF.txt
        open (unit=100,file="POSF.txt",action="write",status="replace")
        write (100,*)  totnum ! , "    ncycle:", ncycle
        write (100,*)  "Factional" 
     ! print final velocity to VELOF.txt
        open (unit=28,file="VELOF.txt",action="write",status="replace")
        write (28,*)  "Velocity of atoms" 
        write (28,*)  totnum  !, "     ncycle:", ncycle
        write (28,*)  "  Reduced unit" 

         do i = 1, totnum
            write (100, "(3F15.7)" ) pos(1,i)/lx,  pos(2,i)/ly, pos(3,i)/lz 
            write (28, "(3F15.7)" )  velo(:,i) 
         end do

       close(100)  ! POSF.txt
       close(28)  ! VELOF.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


write(*,*) "conv press" , 1.0d0/force2reducedforce*length2reducedlength**2

call when

        stop
        end program

