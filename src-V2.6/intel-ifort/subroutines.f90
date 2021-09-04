module subroutines
  use omp_lib
  !Contains all subroutines necessaries in the main program

  !Input data to work in the read routines
  IMPLICIT NONE
  integer :: statusoper
  integer :: tmin, zmin
  real :: poszmin, poszmax, postmin, postmax
  CHARACTER(2),ALLOCATABLE,DIMENSION(:)	::	namematerial		! Names of materials namematerial(1:nummedia,1:2)
  character(len=30) :: filematerialprop,filereflective
  !Quantities which depend on the material, energy and sigma
  REAL,ALLOCATABLE,DIMENSION(:,:,:)			::	vmedia			! Velocity in media in nm/fs
  REAL,ALLOCATABLE,DIMENSION(:,:,:)			::	taumedia		! Lifetime in media in fs
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)	        	::	pmedia			! Prob of abs in media
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)	        	::	pvacuum			! Prob of abs in media

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BEGIN: SubRoutines											!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!													!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate everything						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE allocate_input_variable ()   ! CALLED in readmaterialprop because the number maxnumenergy is needed
    USE global, only: v, tau, p, Seff, flux, n,&
         Emax, zmax, tmax, tbound, nummedia,lftravel,pscatt,Seff_vac_rout
    !Probably we can reduce the dimension in v, tau and p switching tmax from tbound

    ALLOCATE(vmedia(nummedia,2,Emax),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating vmedia'
    ALLOCATE(taumedia(nummedia,2,Emax),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating taumedia'
    ALLOCATE(pmedia(nummedia,2,Emax,2,Emax),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating pmedia'
    ALLOCATE(pvacuum(nummedia,2,Emax,2),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating pvacuum'

    ALLOCATE(v(2,Emax,zmax,tmax*10),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating v'

    ALLOCATE(tau(2,Emax,zmax,tmax),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating tau'

    ALLOCATE(n(2,Emax,zmax,tbound),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating n'

    ALLOCATE(p(2,Emax+1,2,Emax+1,zmax,tmax),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating p'

    ALLOCATE(Seff(2,Emax+1,zmax,tbound),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating Sext'

    ALLOCATE(Seff_vac_rout(zmax,tbound,2),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating Sext_vac_rout'

    ALLOCATE(flux(2,Emax,zmax,tbound),STAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Allocating flux'

    allocate(lftravel(zmax+1,zmax,tmax,zmax), pscatt(zmax+1,zmax,tmax,zmax))
    lftravel(:,:,:,:)=100000000.0
    pscatt(:,:,:,:)=100000000.0
    
  END SUBROUTINE allocate_input_variable



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! General Input							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE readinputgeneral ()
    USE global, only: smooth, k, dz, dt, dE, zmax, tmax, tbound, &
         nummedia, placeinterfaces, extsourcetype, verbosity
    implicit none
    real :: rtbound
    integer :: j,i
    !logical, dimension(:), allocatable :: valid
		 
    WRITE (*,*) '**** SUBROUTINE readinputgeneral'
    WRITE (*,*)
    
    READ (*,*)
    READ (*,*)
    READ (*,*)
    
    READ (*,*) 
    READ (*,*) verbosity
    IF ((verbosity <0).OR.(verbosity>2)) THEN
       WRITE (*,*) 'ERROR: Wrong Verbosity'
    END IF
    
    READ (*,*)
    READ (*,*) k
    
    READ (*,*)
    READ (*,*) smooth
    
    READ (*,*) 
    READ (*,*) poszmin
    
    READ (*,*) 
    READ (*,*) poszmax
    IF (poszmin >= poszmax) THEN
       WRITE (*,*) 'ERROR: zmin >= zmax'
    END IF
    
    READ (*,*) 
    READ (*,*) postmin
    
    READ (*,*) 
    READ (*,*) postmax
    IF (postmin >= postmax) THEN
       WRITE (*,*) 'ERROR: tmin >= tmax'
    END IF    
    
    READ (*,*) 
    READ (*,*) dz
    IF (dz >= (poszmax-poszmin)) THEN
       WRITE (*,*) 'ERROR: dz >= (zmax-zmin)'
    END IF
    !Since the surfaces at borders are dz/2 outside of the material in the calculation. With this shift they will be in the borders 
    poszmax=poszmax-dz  
    
    READ (*,*) 
    READ (*,*) dt
    IF (dt >= (postmax-postmin)) THEN
       WRITE (*,*) 'ERROR: dt >= (tmax-tmin)'
    END IF

    READ (*,*) !For future events, we need take into account events which happened in the past as much as rtbound seconds ago
    READ (*,*) rtbound
    !Number of iteration which are important for future events
    tbound=int(rtbound/dt)
    
    READ (*,*) 
    READ (*,*) dE

    ! Definition of mesh
    zmin=1               !Z_Real = (zpos-1)*dz+dz/2
    zmax=poszmax/dz+1
    tmin=1
    tmax=postmax/dt+1
    
    READ (*,*) 
    READ (*,*) nummedia
    
    ALLOCATE(namematerial(nummedia),STAT=statusoper)
    ALLOCATE(placeinterfaces(0:nummedia),STAT=statusoper)
    placeinterfaces(0)=(zmin-1)*dz
    placeinterfaces(nummedia)=zmax*dz
    
    READ (*,*)
    READ (*,*, IOSTAT=statusoper) placeinterfaces(1:nummedia-1)

    READ (*,*)
    READ (*,*, IOSTAT=statusoper) namematerial
    
    READ (*,*) 
    READ (*,*) filematerialprop,filereflective
    
    READ (*,*) 
    READ (*,*) extsourcetype
    IF ((extsourcetype /= 'step').AND.(extsourcetype /= 'general').AND.(extsourcetype /= 'expgauss')&
         .AND.(extsourcetype /= 'delta')) THEN
       WRITE (*,*) 'ERROR: Wrong Type of input for the External Source'
    END IF
    
    IF ((verbosity ==2).OR.(verbosity ==1)) THEN
       WRITE (*,*) 'verbosity        = ',verbosity
       WRITE (*,*) 'poszmin          = ',poszmin
       WRITE (*,*) 'poszmax          = ',poszmax+1
       WRITE (*,*) 'postmin          = ',postmin
       WRITE (*,*) 'postmax          = ',postmax
       WRITE (*,*) 'dz               = ',dz 
       WRITE (*,*) 'dt               = ',dt
       WRITE (*,*) 'zmin             = ',zmin
       WRITE (*,*) 'zmax             = ',zmax
       WRITE (*,*) 'tmin             = ',tmin
       WRITE (*,*) 'tmax             = ',tmax
       
       WRITE (*,*) 'filematerialprop = ',filematerialprop
       WRITE (*,*) 'filereflective = ',filereflective
       WRITE (*,*) 'extsourcetype    = ',extsourcetype
       WRITE (*,*)
    END IF
    
    ! OUTPUT *********************************************************
    WRITE (*,*) '**** END SUBROUTINE readinputgeneral'		!*
    WRITE (*,*)							!*
    WRITE (*,*)							!*
    ! END OUTPUT *****************************************************
    
  END SUBROUTINE readinputgeneral
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Material Input						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE readmaterialprop ()
    USE global, only: nummedia, Emax, tmax, v, tau, p, dz, zmax, &
         lenghtmatname, placeinterfaces, dE, reflect, verbosity,t0source
    
    INTEGER ::	z,z0,i,t,j,l,z1,z2,t_limit
    integer :: sigma, energy, numtot
    real, dimension(nummedia) :: dupmedia, ddownmedia
    ! From file
    CHARACTER(LEN=lenghtmatname)  ::	matnameread		! Current read material name from database
    REAL  ::	trash		! Used tu put unwanted data
    real, dimension(:,:,:), allocatable :: Rrefpos_left, Rrefpos_right
    logical, dimension(:,:), allocatable :: logic
    logical :: trash_ref
    
    WRITE (*,*) '**** SUBROUTINE readmaterialprop'
    WRITE (*,*)
    
    OPEN (UNIT = 11, FILE = filematerialprop, STATUS = "OLD", ACTION = "READ", IOSTAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Opening filematerialprop'
    OPEN (10, file=filereflective, status = "unknown", IOSTAT=statusoper)
    IF (statusoper /= 0) WRITE (*,*) '! ERROR: Opening filereflective'

    DO i=1, 5, 1
       READ (11,*)						! Simply skips the first 5 lines
    END DO
    READ (11,*) Emax
    
    WRITE (*,*) 'Emin = ',1
    WRITE (*,*) 'Emax = ',Emax
        
    REWIND (UNIT=11)
    CALL allocate_input_variable ()

    DO i=1, nummedia, 1
       DO j=1, 8, 1
          READ (11,*)				! Simply skips the first 8 lines
       END DO
       READ (11,*) matnameread
       DO WHILE ( matnameread /= namematerial(i) )		! WARNING: If the name is not found infinite loop
          DO j=1, (24+6*(Emax-1+2)), 1
             READ (11,*)				! Simply skips the current entry
          END DO
          READ (11,*) matnameread
       END DO
       READ(11,*)
       READ(11,*)
       READ(11,*) dupmedia(i), ddownmedia(i)
       DO j=1, 4, 1					! Simply skips the first 4 lines
          READ (11,*)				! Skips values at 0 energy as well
       END DO
       DO j=1, Emax, 1				! Spin up
          READ (11,*) trash, vmedia(i,2,j), taumedia(i,2,j)
       END DO
       
       DO j=1, 3, 1					! Simply skips the first 4 lines
          READ (11,*)				! Skips values at 0 energy as well
       END DO
       DO j=1, Emax, 1				! Spin down
          READ (11,*) trash, vmedia(i,1,j), taumedia(i,1,j)
       END DO
       READ (11,*)
       DO t=1, 4, 1						! not used in sense of t but to scan spin->spin prob
          DO z0=1, 4, 1					! Simply skips the first 4 lines
             READ (11,*)				! Skips values at 0 energy as well
          END DO
          !If I want to add electron from 0 energy level to the others I have to do this
          !Remove on of the last read(11,*) and add the next one:
          !read(11,*) trash, pvacuum00,pvacuum0(i,MOD(4-t,2)+1,1:Emax,(4-t)/2+1)
          DO z=1, Emax, 1				! not used in sense of z but to scan different energies "from"
             READ (11,*) trash, pvacuum(i,MOD(4-t,2)+1,z,(4-t)/2+1), pmedia(i,MOD(4-t,2)+1,z,(4-t)/2+1,1:Emax)
             pmedia(i,MOD(4-t,2)+1,z,(4-t)/2+1,1:Emax)=pmedia(i,MOD(4-t,2)+1,z,(4-t)/2+1,1:Emax)/dE
             pvacuum(i,MOD(4-t,2)+1,z,(4-t)/2+1)=pvacuum(i,MOD(4-t,2)+1,z,(4-t)/2+1)/dE
          END DO
       END DO
       REWIND (UNIT=11)
    END DO

    numtot=nummedia+1
    !we allocate the reflective vectors
    allocate(reflect%ref_left(2,Emax,zmax+1),reflect%ref_right(2,Emax,zmax+1))
    allocate(reflect%z_sn(zmax+1),reflect%z_dx(zmax+1),reflect%valid(zmax+1))
    ALLOCATE(Rrefpos_left(2,Emax,numtot),Rrefpos_right(2,Emax,numtot))
    allocate(logic(2,nummedia-1))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! WE DEFINE THE REFLECTIVE PARAMETERS FOR ALL THE INTERFACES
    ! WE ALSO DEFINE THE DEXTRA AND SINESTRA (RIGHT AND LEFT) INTERFACES FOR EVERY
    ! POINT IN INTERFACES GRID
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    READ(10,*) !Read the reflective parameters for different energies and spin
    do sigma=1,2
       READ(10,*)  ! Comments SPIN
       READ(10,*)  ! Comments Column meaning
       read(10,*) trash_ref,logic(sigma,1:nummedia-1),trash_ref
       do energy=1,Emax
          READ (10,*) trash,Rrefpos_left(sigma,energy,1:numtot)
          READ (10,*) trash,Rrefpos_right(sigma,energy,1:numtot)
       end do
    end do
    reflect%ref_left(:,:,:)=0.0
    reflect%ref_right(:,:,:)=0.0
    reflect%valid(:) = .FALSE.
    ! left end
    reflect%ref_left(:,:,1)=Rrefpos_left(:,:,1)  
    reflect%ref_right(:,:,1)=Rrefpos_right(:,:,1)
    reflect%valid(1)= .TRUE.  ! originally, this was set to TRUE
    ! right end
    reflect%ref_left(:,:,zmax+1)=Rrefpos_left(:,:,nummedia+1)
    reflect%ref_right(:,:,zmax+1)=Rrefpos_right(:,:,nummedia+1)
    reflect%valid(zmax+1)= .TRUE.  ! originally, this was set to TRUE
    !Asignation of the reflective parameters at the spatial grid, i.e.
    !each reflective parameter is assigned to an spatial position
    do j=1,nummedia-1
       if(logic(1,j).and.logic(2,j))then
          i=int(placeinterfaces(j)/dz)+1    !Integer number for the surface
          reflect%ref_left(:,:,i)=Rrefpos_left(:,:,j+1)
          reflect%ref_right(:,:,i)=Rrefpos_right(:,:,j+1)
          reflect%valid(i)= .TRUE.
       end if
    end do

    do sigma=1,2
       do energy=1,Emax
          do i=1,zmax+1
             reflect%z_sn(i)=i
             reflect%z_dx(i)=i
             do j=i-1,1,-1
                reflect%z_sn(i)=j
                if(reflect%valid(j))exit
             end do
             do j=i+1,zmax+1
                reflect%z_dx(i)=j
                if(reflect%valid(j)) exit
             end do
          end do  !Close loop on positions
       end do     !Close loop on energies
    end do        !Close loop on spins

    !If to check the reflective parameters
!    if(verbosity==2)then
!       do sigma=1,2
!          do energy=1,Emax
!             do i=1,zmax+1
!                write(*,*)sigma,energy,i,reflect%z_sn(i),reflect%ref_left(sigma,energy,reflect%z_sn(i)),reflect%z_dx(i),reflect%ref_right(sigma,energy,reflect%z_dx(i))
!             end do
!          end do
!       end do
!    end if
    DEALLOCATE(Rrefpos_left,Rrefpos_right)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CLOSE(11)
    close(10)
    
    !dpop=0

    z0=zmin-1					! Not used in the sense of z0
    i=1						! Inizialization For-loop on i
    DO WHILE (i<=nummedia)						! For-loop on i
       IF (i==nummedia) THEN
          z=zmax
       ELSE
          z=INT(placeinterfaces(i)/dz)
       END IF
       DO j=1, 2, 1					! Loop on spin channels
          DO t=1, Emax, 1	! Loop on energies "from"
             t_limit=t0source+int(taumedia(i,j,t))
             if(int(taumedia(i,j,t))>70) t_limit=t0source+70
             v(j,t,(z0+1):z,1:t_limit)=vmedia(i,j,t)/1.0
             v(j,t,(z0+1):z,t_limit+1:tmax)=vmedia(i,j,t)/1.0
             tau(j,t,(z0+1):z,1:tmax)=taumedia(i,j,t)
             DO l= 1, Emax, 1	! Loop on energy "to"
                p(1,t,j,l,(z0+1):z,1:tmax)=pmedia(i,j,t,1,l)
                p(2,t,j,l,(z0+1):z,1:tmax)=pmedia(i,j,t,2,l)
                if(l.eq.1)then
                   p(1,t,j,Emax+1,(z0+1):z,1:tmax)=pvacuum(i,j,t,1)
                   p(2,t,j,Emax+1,(z0+1):z,1:tmax)=pvacuum(i,j,t,2)
                end if
             END DO
          END DO
       END DO
       !M(tmin-1,(z0+1):z)=dupmedia(i)-ddownmedia(i)
       !dpop(1,tmin-1,(z0+1):z)=dupmedia(i)
       !dpop(0,tmin-1,(z0+1):z)=ddownmedia(i)
       z0=z
       i=i+1
    END DO
    DEALLOCATE(vmedia,taumedia,pmedia,pvacuum)
    
    WRITE (*,*) '**** END SUBROUTINE readmaterialprop'
    WRITE (*,*)
    WRITE (*,*)
  END SUBROUTINE readmaterialprop
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! External Source Input						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE readextsource ()
    USE global, only: extsourcetype, lambdacoefficients, excavrel, &
         timewidth, t0source, Emax, nummedia
    
    INTEGER					::	z,t,i,z1,z2,j,l
    CHARACTER(LEN=3)			::	energystr
    
    WRITE (*,*) '**** SUBROUTINE readextsource'
    WRITE (*,*)
    
    IF (extsourcetype=='delta') THEN
       WRITE (*,*) 'Computing extsourcetype=delta'
    ELSE IF (extsourcetype=='test') THEN
       WRITE (*,*) 'Computing extsourcetype=test'
    ELSE IF (extsourcetype=='expgauss') THEN
       WRITE (*,*) 'Computing expgauss'
       ALLOCATE(lambdacoefficients(nummedia),STAT=statusoper)
       ALLOCATE(excavrel(2,Emax,nummedia),STAT=statusoper)
       DO i=1, 3, 1
          READ (*,*)
       END DO
       READ (*,*) timewidth
       READ (*,*)
       READ (*,*) t0source
       READ (*,*)
       READ (*,*) (lambdacoefficients(i),i=1,nummedia)
       READ (*,*)
       READ (*,*)
       DO i=1, Emax, 1
          READ (*,*) excavrel(1,i,1:nummedia)
       END DO
       READ (*,*)
       DO i=1, Emax, 1
          READ (*,*) excavrel(2,i,1:nummedia)
       END DO
    ELSE IF (extsourcetype=='general') THEN
       !READ (11,*, IOSTAT=statusoper) Sext
    END IF
     
    WRITE (*,*) '**** END SUBROUTINE readextsource'
    WRITE (*,*)
    WRITE (*,*)
  END SUBROUTINE readextsource
  


  ! Routine to calculate Seff(sigma,enrgy,z',t')
  !     Seff(sigma,energy,z',t')=Sscat(sigma,enrgy,z',t')n(sigma,enrgy,z',t')+Sexternal(sigma,enrgy,z',t'+1)
  !         Sscat(sigma,enrgy,z',t')n(sigma,energy,z',t')=sum_sigmaprim int_energyptim p()/tau()*n()
  !         Sexternal(sigma,enrgy,z',t') is a function (the external source)
  ! 
  !   IMPORTANT: THE EXTERNAL SOURCE IS EVALUATED IN TIME=TIME' + 1.
  !-------------------------------------------------------------------------------------------------------------
  subroutine Seffevaluat(sigma,energy,zposprim,timeprim)
    use global, only: Emax, n, tau, p, Seff, prop, dE, dt,extsourcetype
    use functions, only: Sexternal, rectangular, funcmod
    implicit none
    real Ieprim, Sext, IEninteg
    real, dimension(:), allocatable :: QEninteger
    integer :: energy, sigma, enerprim, sigprim,tscale
    integer, intent(in) :: zposprim, timeprim

    allocate(QEninteger(Emax))
    tscale=funcmod(timeprim)
    do enerprim=1,Emax
       QEninteger(enerprim)=0.0
       do sigprim=1,2
            QEninteger(enerprim)=QEninteger(enerprim)+p(sigprim,enerprim,sigma,energy,zposprim,timeprim)&              
                 *n(sigprim,enerprim,zposprim,tscale)*(1.0-exp(-dt/tau(sigprim,enerprim,zposprim,timeprim)))
       end do
    end do
    IEninteg=rectangular(dE,Emax,QEninteger)  !Rectangular method to compute the integral
    Sext=Sexternal(sigma,energy,zposprim,timeprim+1)
    Ieprim=IEninteg+Sext
    Seff(sigma,energy,zposprim,tscale)=Ieprim
    if(energy==(Emax+1)) Seff(sigma,energy,zposprim,tscale)=IEninteg
    deallocate(QEninteger)  
    return
  end subroutine Seffevaluat



  !Function to compute the electron probability of being scattered and to go to the Marco's sea
  !This is computed for every energy, spin and for every position.
  subroutine Seffvacuum_rut(zposprim,timeprim,sigma)
    use global, only: Emax, n, tau, p, dE,dt,Seff_vac_rout
    use functions, only: rectangular, funcmod
    implicit none
    real IEninteg
    real, dimension(:), allocatable :: QEninteger
    integer :: enerprim, sigprim, tscale
    integer, intent(in) :: zposprim, timeprim,sigma
         
    allocate(QEninteger(Emax))
    tscale=funcmod(timeprim)
    do enerprim=1,Emax
       QEninteger(enerprim)=0.0
       do sigprim=1,2
          QEninteger(enerprim)=QEninteger(enerprim)+p(sigprim,enerprim,sigma,Emax+1,zposprim,timeprim)&
               *n(sigprim,enerprim,zposprim,tscale)*(1.0-exp(-dt/tau(sigprim,enerprim,zposprim,timeprim)))          
       end do
    end do
    IEninteg=rectangular(dE,Emax,QEninteger)
    Seff_vac_rout(zposprim,tscale,sigma)=IEninteg
    deallocate(QEninteger)
    
    return
  end subroutine Seffvacuum_rut



  !Routine which evaluates the electron's lifetime without scattering in a given 
  !position (z) and time (t), and by depending on the z_0 value. If the 
  !electron velocity is constant this quantity is also constant, 
  !but if v=v(z_0) we have to do an explicit integral
  !         lftravel(z,t,z_0)=int_{z_0}^{z} dz'/v(z')
  !
  !It also evaluates the probability of scattering of an electron in a given
  !position (z) and time (t), and by depending on the z_0 value. If the 
  !lifetive and electron velocity are constants this quantity is also constant, 
  !but if v=v(z_0) or tau=tau(z_0) we have to do an explicit integral
  !         pscatt(z,t,z_0)=int_{z_0}^{z} dz'/(v(z')*tau(z'))
  ! 
  !By taking into account the sccatering both the life time and the probability of 
  !scattering are evaluated as sums as follow:
  !
  !   lftravel(z,t,z_0)=lftravel(z,t,z_0+1)+dz/v(Sigma,E,z_0,t_0)
  !with
  !   t_0=t-int(lftravel(z_0+1)/dt)-1
  !
  !i.e. we evaluate both quantities as the sumation of these quantities in the 
  !travel from one point to the next one.
     !check boundaries

  Subroutine timetravel(z_half,R_t)
    use global, only: dz, dt, iint, v, tau, zmax, verbosity, &
         reflect,last_z0_pforz,last_z0_mforz,lftravel,pscatt
    implicit none
    integer, intent(in) :: z_half
    integer :: zprime, tprime, new_time,z_pos
    integer, pointer :: zdx_point,zsn_point
    real, intent(in) :: R_t    

    new_time=int(R_t/dt)+1   !to avoid new_time==0 we add +1
    !_____________________________________________________________________________________________________
    !Calculation of the time and probability that an electron takes to go from the interface "z_half" to 
    !the point "zpos" without scattering. In this case is calculated from left to right(positive quantities)
    !_____________________________________________________________________________________________________
    z_pos=z_half-1
    if(z_half > 1)then
       lftravel(z_half,z_pos,new_time,iint%zpos)=dz/(2.0*v(iint%sigma,iint%energy,z_pos,new_time))
       pscatt(z_half,z_pos,new_time,iint%zpos)=dz/(2.0*v(iint%sigma,iint%energy,z_pos,new_time)*&
            tau(iint%sigma,iint%energy,z_pos,new_time))
       tprime=int((R_t-lftravel(z_half,z_pos,new_time,iint%zpos))/dt)+1
    else
       tprime=1
    end if
    !Calculation further into the left point from which the electron can come
    zsn_point=>reflect%z_sn(z_half)
    if(tprime.lt.2) then
       last_z0_mforz=zmax+1 !z_pos+1  
    else if(tprime >= 2.and.(zsn_point > (z_pos-1)))then 
       last_z0_mforz=z_pos
    else
       do zprime=z_pos-1,zsn_point,-1
          lftravel(z_half,zprime,new_time,iint%zpos)=lftravel(z_half,zprime+1,new_time,iint%zpos)+dz/(2*v(iint%sigma,iint%energy,zprime+1,tprime))+&
               dz/(2*v(iint%sigma,iint%energy,zprime,tprime))
          pscatt(z_half,zprime,new_time,iint%zpos)=pscatt(z_half,zprime+1,new_time,iint%zpos)+dz/(2*v(iint%sigma,iint%energy,zprime+1,tprime)*&
               tau(iint%sigma,iint%energy,zprime+1,tprime))+dz/(2*v(iint%sigma,iint%energy,zprime,tprime)*&
               tau(iint%sigma,iint%energy,zprime,tprime))
          tprime=int((R_t-lftravel(z_half,zprime,new_time,iint%zpos))/dt)+1
          if(tprime.lt.2) then
             last_z0_mforz=zprime+1
             exit
          end if
          last_z0_mforz=zsn_point
       end do
    end if
    
    !_____________________________________________________________________________________________________
    !Calculation of the time and probability that an electron takes to go from the interface "z_half" to 
    !the point "zpos" without scattering. In this case is calculated from right to left(negative quantities)
    !_____________________________________________________________________________________________________    
    z_pos=z_half
    if(z_half <= zmax)then
       lftravel(z_half,z_pos,new_time,iint%zpos)=-dz/(2*v(iint%sigma,iint%energy,z_pos,new_time))
       pscatt(z_half,z_pos,new_time,iint%zpos)=-dz/(2*v(iint%sigma,iint%energy,z_pos,new_time)*tau(iint%sigma,iint%energy,z_pos,new_time))
       tprime=int((R_t+lftravel(z_half,z_pos,new_time,iint%zpos))/dt)+1
    else
       tprime=1
    end if
    !Calculation further into the right point from which the electron can come
    zdx_point=>reflect%z_dx(z_half)
    if(tprime.lt.2) then
       last_z0_pforz=-1 !z_pos-1   
    else if(tprime >= 2.and.((zdx_point-1) < (z_pos+1)))then
       last_z0_pforz=z_pos
    else
       do zprime=z_pos+1,zdx_point-1,+1
          lftravel(z_half,zprime,new_time,iint%zpos)=lftravel(z_half,zprime-1,new_time,iint%zpos)-&
               dz/(2*v(iint%sigma,iint%energy,zprime-1,tprime))-dz/(2*v(iint%sigma,iint%energy,zprime,tprime))
          pscatt(z_half,zprime,new_time,iint%zpos)=pscatt(z_half,zprime-1,new_time,iint%zpos)-&
               dz/(2*v(iint%sigma,iint%energy,zprime-1,tprime)*tau(iint%sigma,iint%energy,zprime-1,tprime))&
               -dz/(2*v(iint%sigma,iint%energy,zprime,tprime)*tau(iint%sigma,iint%energy,zprime,tprime))
          tprime=int((R_t+lftravel(z_half,zprime,new_time,iint%zpos))/dt)+1
          if(tprime.lt.2) then
             last_z0_pforz=zprime-1
             exit
          end if
          last_z0_pforz=zdx_point-1
       end do
    end if

  end Subroutine timetravel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine phi_right(repet_index,z_half,t_iter,delta_t,delta_t_tau,value_phi)
    use global, only: reflect, iint, dt, dz, v, tau, tbound, Seff,Rtau,&
         last_z0_pforz,last_z0_mforz,lftravel,pscatt,Seff_vac_rout
    use functions, only: Phi,funcmod
    real :: life_time,prob_scatt,R_t,t_m_delta
    real :: delta_t_new, delta_t_tau_new,value_phi_new
    real :: Atmp
    real, intent(in) :: delta_t, delta_t_tau
    real, intent(inout) :: value_phi
    integer :: t, tscale, tprime, zprime,z_pos_z0,repet_index_new
    integer, intent(in) :: t_iter,repet_index
    integer :: z_half,t_low_bound
    integer, pointer :: p_left

    R_t=(t_iter-1)*dt                                               !We calculate the real time
    t_m_delta=R_t-delta_t   !! maximum time?? 
    if(delta_t < 0.0) then
       write(*,*) 'ERROR:  delta_t < 0 for particles coming from the right side' 
       stop
    end if   

    if(abs(delta_t).lt.(tbound-1)*dt.and.t_m_delta.gt.0.0)then !! I do not understand this condition
       !_____________________________________________________________________________________________________________________________
       !Calculation of the lifetime and scattering probability of electrons from the interface "z_half" to all the spatial points
       !_____________________________________________________________________________________________________________________________
       call timetravel(z_half,t_m_delta)
       p_left=>reflect%z_sn(z_half)                                                        !Pointer to  reflect%z_sn(z_half) 
       !$OMP PARALLEL IF(t_iter > 20) PRIVATE(t,t_low_bound,tscale,Atmp) !!!REDUCTION(+:value_phi)
       !----------------------------------------
       ! calculate left-hand side of Eq. (18)
       !----------------------------------------
       Atmp=0.0
       !$OMP DO 
!       do z_pos_z0=p_left,z_half-1
!          if((last_z0_mforz > z_pos_z0).or.(z_pos_z0 > last_z0_pforz))cycle               !We avoid points which don't have time to reach to z
       do z_pos_z0=last_z0_mforz,z_half-1
          t_low_bound=t_iter-tbound
          if(t_low_bound <= 1)t_low_bound=2
          do t=t_low_bound,t_iter
             tscale=funcmod(t-1)
             Atmp=Atmp+Seff(iint%sigma,iint%energy,z_pos_z0,tscale)*Phi(z_half,z_pos_z0,t,delta_t,delta_t_tau)  !! S^e * psi
          end do
       end do
       !$OMP END DO
       !$OMP ATOMIC               !We do not need the atomic clause anymore because the sum is done with the reduction clause     
       value_phi=value_phi+Atmp
       !$OMP END PARALLEL
       !
       if(last_z0_mforz <= p_left)then!.and.p_left <= last_z0_pforz)then                 !We avoid points which don't have time to reach to z
          tprime=int(t_m_delta/dt)+1
          life_time=lftravel(z_half,p_left,tprime,iint%zpos)                  !Calculation of life_time
          prob_scatt=pscatt(z_half,p_left,tprime,iint%zpos)                   !Calculation of prob_scatt

          tprime=int((t_m_delta-life_time)/dt)+1
          delta_t_new=delta_t+life_time+dz/(2.0*v(iint%sigma,iint%energy,p_left,tprime))              !Calculation of the new delta_t [Eq.(18)]
          delta_t_tau_new=delta_t_tau+prob_scatt+dz/(2.0*v(iint%sigma,iint%energy,p_left,tprime)*&    !Calculation of the new delta_t_tau [Eq. (18)]
               tau(iint%sigma,iint%energy,p_left,tprime))
       else
          delta_t_new=R_t+1.0
          delta_t_tau_new=0.0
       end if

       if(p_left.ne.1)then
          value_phi_new=0.0
          repet_index_new=2*repet_index+1
          call phi_right(repet_index_new,p_left,t_iter,delta_t_new,delta_t_tau_new,value_phi_new)
          value_phi=value_phi+(1.0-reflect%ref_right(iint%sigma,iint%energy,p_left))*value_phi_new
       end if
       delta_t_new=-delta_t_new
       delta_t_tau_new=-delta_t_tau_new
       value_phi_new=0.0
       repet_index_new=2*repet_index+2
       call phi_left(repet_index_new,p_left,t_iter,delta_t_new,delta_t_tau_new,value_phi_new)
       value_phi=value_phi-reflect%ref_left(iint%sigma,iint%energy,p_left)*value_phi_new
    end if
    !if(value_phi>1E-06.and.repet_index>2000)print*, 'right',repet_index,delta_t,'s',delta_t_tau,last_z0_mforz,last_z0_pforz,value_phi
    
  end subroutine phi_right
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  recursive subroutine phi_left(repet_index,z_half,t_iter,delta_t,delta_t_tau,value_phi)
    use global, only: reflect, iint, dt, dz, v, tau, tbound, Seff,&
         zmax,last_z0_pforz,last_z0_mforz, lftravel, pscatt,Seff_vac_rout
    use functions, only: Phi,funcmod
    real :: life_time,prob_scatt,R_t, t_m_delta
    real :: delta_t_new, delta_t_tau_new,value_phi_new
    real :: Atmp
    real, intent(in) :: delta_t, delta_t_tau
    real, intent(inout) :: value_phi
    integer :: t, tscale,tprime, zprime, z_pos_z0,repet_index_new
    integer, intent(in) :: t_iter,repet_index
    integer :: z_half,t_low_bound
    integer, pointer :: p_right

    R_t=(t_iter-1)*dt                                                          !We calculate the real time
    t_m_delta=R_t+delta_t
    if(delta_t > 0.0) then
       write(*,*) 'ERROR:  delta_t > 0 for particles coming from the left side' 
       stop
    end if

    if(abs(delta_t).lt.(tbound-1)*dt.and.t_m_delta.gt.0.0)then
       !_____________________________________________________________________________________________________________________________
       !Calculation of the lifetime and scattering probability of electrons from the interface "z_half" to all the the spatial points
       !_____________________________________________________________________________________________________________________________
       call timetravel(z_half,t_m_delta)
       p_right=>reflect%z_dx(z_half)                                 !Pointer to  reflect%z_sn(z_half) 
       !$OMP PARALLEL IF(t_iter > 20) PRIVATE(t,t_low_bound,tscale,Atmp) !!!REDUCTION(+:value_phi)
       Atmp=0.0
       !$OMP DO 
!       do z_pos_z0=z_half,p_right-1
!          if(last_z0_mforz > z_pos_z0.or.z_pos_z0 > last_z0_pforz)cycle                 !We avoid points which don't have time to reach to z
       do z_pos_z0=z_half,last_z0_pforz
          t_low_bound=t_iter-tbound
          if(t_low_bound <= 1)t_low_bound=2
          do t=t_low_bound,t_iter
             tscale=funcmod(t-1)
             Atmp=Atmp+Seff(iint%sigma,iint%energy,z_pos_z0,tscale)*Phi(z_half,z_pos_z0,t,delta_t,delta_t_tau)
          end do
       end do
       !$OMP END DO
       !$OMP ATOMIC               !We do not need the atomic clause anymore because the sum is done with the reduction clause     
       value_phi=value_phi+Atmp
       !$OMP END PARALLEL

       if((p_right-1).le.last_z0_pforz)then !last_z0_mforz.le.(p_right-1).and.                !We avoid points which don't have time to reach to z
          tprime=int(t_m_delta/dt)+1
          life_time=lftravel(z_half,p_right-1,tprime,iint%zpos)                  !Calculation of life_time
          prob_scatt=pscatt(z_half,p_right-1,tprime,iint%zpos)                   !Calculation of prob_scatt
          
          tprime=int((t_m_delta+life_time)/dt)+1
          delta_t_new=delta_t+life_time-dz/(2.0*v(iint%sigma,iint%energy,p_right-1,tprime))              !Calculation of the new delta_t
          delta_t_tau_new=delta_t_tau+prob_scatt-dz/(2.0*v(iint%sigma,iint%energy,p_right-1,tprime)*&    !Calculation of the new delta_t_tau
               tau(iint%sigma,iint%energy,p_right-1,tprime))
       else
          delta_t_new=-R_t-1.0
          delta_t_tau_new=0.0
       end if

       if(p_right < zmax+1)then
          value_phi_new=0.0
          repet_index_new=2*repet_index+2
          call phi_left(repet_index_new,p_right,t_iter,delta_t_new,delta_t_tau_new,value_phi_new)
          value_phi=value_phi+(1.0-reflect%ref_left(iint%sigma,iint%energy,p_right))*value_phi_new
       end if
       delta_t_new=-delta_t_new
       delta_t_tau_new=-delta_t_tau_new
       value_phi_new=0.0
       repet_index_new=2*repet_index+1
       call phi_right(repet_index_new,p_right,t_iter,delta_t_new,delta_t_tau_new,value_phi_new)
       value_phi=value_phi-reflect%ref_right(iint%sigma,iint%energy,p_right)*value_phi_new 
    end if
    !if(value_phi>1E-05.and.repet_index>2000)print*, 'left',repet_index,delta_t,'s',delta_t_tau,last_z0_mforz,last_z0_pforz,value_phi

  end subroutine phi_left



  Subroutine write_files(succesful)
    use global, only: Emax
    implicit none
    integer :: E_index
    logical :: succesful
    CHARACTER(LEN=3)   :: energystr

    DO E_index=1, Emax
       WRITE ( energystr( 1:3 ), '(I3.3)'), E_index
       OPEN (UNIT = 20+E_index, FILE = &
            './output/ndownE'// energystr(1:3) //'.out', &
            ACTION = "WRITE", IOSTAT=statusoper)
       if(statusoper.ne.0)succesful=.FALSE.
       OPEN (UNIT = 40+E_index, FILE = &
            './output/nupE'// energystr(1:3) //'.out', &
            ACTION = "WRITE", IOSTAT=statusoper)
       if(statusoper.ne.0)succesful=.FALSE.
    END DO
    open(98, FILE = './output/ndownE000.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(99, FILE = './output/nupE000.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(621, FILE = './output/ntotal_down.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(622, FILE = './output/ntotal_up.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(623, FILE = './output/fluxtotal_down.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(624, FILE = './output/fluxtotal_up.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(625, FILE = './magnetization.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(626, FILE = './flux.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(721, FILE = './output/js_right_down.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(722, FILE = './output/js_right_up.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(723, FILE = './output/js_left_down.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.
    open(724, FILE = './output/js_left_up.out',ACTION = "WRITE", IOSTAT=statusoper)
    if(statusoper.ne.0)succesful=.FALSE.

    return
  end Subroutine write_files


  subroutine SC_energy_dependent(valid,sc_right,sc_left)
    use global, only: Emax,iint,zmax,tmax,JS_total
    implicit none
    integer :: E_index
    real, intent(in) :: sc_right,sc_left
    !real, dimension(Emax,zmax,tmax) :: JS_total
    character(len=20)  :: command,directory
    CHARACTER(LEN=3)   :: energystr
    logical :: valid

    if(valid) then                 !Create directory
       allocate(JS_total(Emax,zmax,tmax))
       directory='spin_current'
       command='mkdir '//directory//''
       CALL system(command)
       
       DO E_index=1, Emax
          WRITE ( energystr( 1:3 ), '(I3.3)'), E_index
          OPEN (UNIT = 200+E_index, FILE = &
               './spin_current/js_total_E'// energystr(1:3) //'.dat', &
               ACTION = "WRITE", IOSTAT=statusoper)
       END DO
    else                           !Save data
       if(iint%sigma==1)then
          JS_total(iint%energy,iint%zpos,iint%time)=0.0
          JS_total(iint%energy,iint%zpos,iint%time)=JS_total(iint%energy,iint%zpos,iint%time)-sc_right-sc_left
       else
          JS_total(iint%energy,iint%zpos,iint%time)=JS_total(iint%energy,iint%zpos,iint%time)+sc_right+sc_left
          E_index=iint%energy
          WRITE (200+E_index,'(F18.8,$)') JS_total(iint%energy,iint%zpos,iint%time)
          if(iint%zpos==zmax)WRITE (200+E_index,*) 
       end if
    end if

    return
  end subroutine SC_energy_dependent


  

  !Routine which check the different important values for the program, like 
  !           Seff, Sexternal, Top, v, tau, lifetime, pscatt
  !It works by writing this values in several output files
  Subroutine checkitems(skip)
    use global, only: Seff, v, tau, iint, Emax, zmax,tmax,flux
    use functions, only: Sexternal, funcmod
    implicit none
    integer :: tscale, S_index, E_index, skip
    CHARACTER(3) :: energystr
   
    if(iint%sigma == iint%energy.and.iint%energy == iint%zpos.and.iint%zpos == 1.and.iint%time == 2)then
       DO E_index=1, Emax          
          WRITE ( energystr( 1:3 ), '(I3.3)'), E_index
          OPEN (UNIT = 120+E_index, FILE = './output/taudownE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 140+E_index, FILE = './output/tauupE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 220+E_index, FILE = './output/vdownE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 240+E_index, FILE = './output/vupE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 320+E_index, FILE = './output/SeffdownE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 340+E_index, FILE = './output/SeffupE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 420+E_index, FILE = './output/SextdownE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 440+E_index, FILE = './output/SextupE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 520+E_index, FILE = './output/fluxdownE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
          OPEN (UNIT = 540+E_index, FILE = './output/fluxupE'// energystr(1:3) //'.out',ACTION = "WRITE", IOSTAT=statusoper)
       END DO
    end if

    tscale=funcmod(iint%time-1)
    WRITE (120+(iint%sigma-1)*20+iint%energy,'(F18.8,$)')tau(iint%sigma,iint%energy,iint%zpos,iint%time-1)
    WRITE (220+(iint%sigma-1)*20+iint%energy,'(F18.8,$)')v(iint%sigma,iint%energy,iint%zpos,iint%time-1)
    WRITE (320+(iint%sigma-1)*20+iint%energy,'(F18.8,$)')Seff(iint%sigma,iint%energy,iint%zpos,tscale)
    WRITE (420+(iint%sigma-1)*20+iint%energy,'(F18.8,$)')Sexternal(iint%sigma,iint%energy,iint%zpos,iint%time-1)
    WRITE (520+(iint%sigma-1)*20+iint%energy,'(F18.8,$)')flux(iint%sigma,iint%energy,iint%zpos,tscale)

    !Skip a line after writing for each time every Sigma, Energy and position
    if(skip.eq.1) then
       do S_index=1,2
          do E_index=1,Emax
             WRITE (120+(S_index-1)*20+E_index,*) 
             WRITE (220+(S_index-1)*20+E_index,*) 
             WRITE (320+(S_index-1)*20+E_index,*) 
             WRITE (420+(S_index-1)*20+E_index,*) 
             WRITE (520+(S_index-1)*20+E_index,*) 
          end do
       end do
    end if

!    DO E_index=1, Emax
!       CLOSE (UNIT = 20+E_index)
!       CLOSE (UNIT = 40+E_index)
!    END DO

    return
  end Subroutine checkitems
  
  
  !!													!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! END: Subroutines											!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module subroutines
