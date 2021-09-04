MODULE global
use omp_lib
! Contains all data regarding the solution of the equation

  IMPLICIT NONE
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Input Variable						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Define Variables type: --> vector
        !                       --> prop
        type vector 
           real :: vecval
        end type vector

        type prop
           integer :: sigma
           integer :: energy
           integer :: zpos
           integer :: time
        end type prop
        
        type reflection
           real, dimension(:,:,:), allocatable :: ref_left, ref_right
           integer, dimension(:), allocatable :: z_sn,z_dx
           logical, dimension(:), allocatable :: valid
        end type reflection
        
        type (prop), target :: iint
        type (reflection), target :: reflect

        integer :: last_z0_pforz,last_z0_mforz

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Solution Variable
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL, DIMENSION(:,:,:,:), allocatable	::	n		! 1/nm		! array totaln[2,Emax,zmax,tmax]

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Singletime Variables						!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! STORED
        real, dimension(:,:,:,:), allocatable :: Seff   ! 1/(nm*fs)       stored array Seff[2,Emax,zmax,tmax]
        real, dimension(:,:,:,:), allocatable :: flux   ! 1/(nm*fs)       stored array flux[2,Emax,zmax,tmax]
        !type (vector), dimension(:,:,:,:), allocatable :: Qen    ! 1/(nm*fs*eV)    stored array Qen[2,Emax,zmax,tmax]
        real, dimension(:,:,:), allocatable :: Seff_vac_rout   ! 1/(nm*fs)       stored array Seff[2,Emax,zmax,tmax]
        real, dimension(:,:,:), allocatable :: JS_total ! 1/(nm*fs)        stored array JS_total[Emax,zmax,tmax]

        
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Variables for save the deltat and deltattau for each S,E,z,t	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        real, dimension(:,:,:,:), allocatable :: lftravel, pscatt


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Input Variables: They are readed in subroutines read_...	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer :: verbosity
		INTEGER,PARAMETER					::	lenghtmatname=2
		CHARACTER(LEN=10)					::	extsourcetype		
        !		  'step' if fileextsource is in step format	   												      ! 		'general' if it is in numerical format
		!		  'expgauss' if it is gaussian in time and exp in space
		!		  'test' for test source

        !Mesh for position, time and Energy. Variables which allocate the arrays
		REAL 			                          	::	dz              ! nm		! Step Spatial Mesh
		REAL                         				::	dt              ! fs		! Step Temporal Mesh
        REAL                                        ::  dE              ! eV            ! Step Energy Mesh
		INTEGER                         			::	zmax	
		INTEGER                         			::	tmax	
        INTEGER                                     ::  tbound
		INTEGER                        				::	Emax

        !Materials
		INTEGER							::	nummedia		! Number of media
        REAL			                ::	k		        ! nm		! smearing theta function
       	REAL			                ::	smooth

        !Quantities which depend on the material
		REAL,ALLOCATABLE,DIMENSION(:)				::	placeinterfaces		! Position of Interfaces in nm
		REAL,ALLOCATABLE,DIMENSION(:)				::	lambdacoefficients	! lambdacoefficients[1:nummedia]			

                !Quantities which depend on the material, energy and sigma
		REAL,ALLOCATABLE,DIMENSION(:,:,:)			::	excavrel		! exc_average_el[2,Emax,nummedia]

                !Variables related with the external source (extysourcetype). excavrel is also related with the external source
		REAL							::	timewidth
		REAL							::	t0source		! fs		! time position peak source

        	! Material properties (will be updated in the run!!!)
		REAL,ALLOCATABLE,DIMENSION(:,:,:,:)	                ::	v		! nm/fs		! array v[2,Emax,zmax,tmax]	
        REAL, DIMENSION(:,:,:,:), allocatable	            ::	tau		! fs		! array tau[2,Emax,zmax,tmax]
		REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:,:)	            ::	p		! Prob transf	! array p[2,Emax,2,Emax,zmax,tmax]

        REAL                                                ::  Rtau    ! tau for a given Sigma,E,z and time

END MODULE global

