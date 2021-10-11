  program interface
    use omp_lib
    use global, only: tmax, zmax, Emax, iint, reflect, n, tau, v,&
         Seff, flux, dt, dz, tbound, verbosity,Rtau !,pscatt, lftravel 
    use subroutines
    use functions
    implicit none

    real*8 :: verystarttime,starttime,endtime
    integer :: TID, MAX, number_of_threads
    integer :: t_index,z_index,E_index,S_index                      !Indixes for time, position, energy and spin
    integer :: timeprim, zposprim                                   !Indixes for alternative time and position
    integer :: timem1,tscale,tscale1                                !Time variables  timem1=t_index-1;  tscale=mod(timeprim,tbound)
    integer, pointer :: z_ref,p_left
    real :: Item, MS                                                !Real variables which for working in the // operator
    real :: delta_t,delta_t_tau
    real :: value_phi,flux_value                                    !Real variables to keep n()%vecval and tau()%vecval
    real :: value_phi_next_right,value_phi_next_left
    real :: sc_right,sc_left
    real, dimension(:,:,:), allocatable :: nvacuum                  !Number of electrons in the Marco's sea
    real, dimension(:,:,:), allocatable :: n_total                  !Number of total electrons with spin up and down
    real, dimension(:,:,:), allocatable :: spin_current_right,spin_current_left
    real, dimension(:,:,:), allocatable :: flux_total               !Total flux  electrons with spin up and down
    real, dimension(:,:), allocatable :: magnetization              !Magnetization for each t and z (it is the difference between e- up and down including Marco's see)
    real, dimension(:,:), allocatable :: flux_F                     !Total flux of electrons
    character(len=20)  :: command,directory
    logical :: succesful

    directory='output'
    command='mkdir '//directory//''
    succesful=.TRUE.

    !It tells to the program the number of threads ("processors") to use 
    !in the parallel version
    !It sets the number of threads to be used by subsequent parallel regions
    !number_of_threads=2
    !call omp_set_num_threads(number_of_threads)

    ! Get Input
    CALL readinputgeneral ()
    CALL readmaterialprop ()
    CALL readextsource ()    

    !Allocata some vector which will only be used in the main rutine
    !allocate(lftravel(zmax+1,zmax,tmax), pscatt(zmax+1,zmax,tmax))
    allocate(nvacuum(2,zmax,tmax),n_total(2,zmax,tmax),flux_total(2,zmax,tmax))
    allocate(spin_current_right(2,zmax,tmax),spin_current_left(2,zmax,tmax))
    allocate(magnetization(zmax,tmax),flux_F(zmax,tmax))

    !Initialite the number of electron and Seff vectors
    Seff(1:2,1:Emax,1:zmax,1)=0.0
    n(1:2,1:Emax,1:zmax,1:tbound)=0.0

    !Open the files where will be written the number of electrons for each energy and spin.
    CALL system(command)
    call write_files(succesful)

    nvacuum(1:2,1:zmax,1)=0.0
    n_total(1:2,1:zmax,1)=0.0
    spin_current_left(1:2,1:zmax,1:tmax)=0.0
    spin_current_right(1:2,1:zmax,1:tmax)=0.0
    flux_total(1:2,1:zmax,1)=0.0
    magnetization(:,:)=0.0
    flux_F(:,:)=0.0
    verystarttime = omp_get_wtime()
    call SC_energy_dependent(.TRUE.,sc_right,sc_left)

    !Start the main program with loops over time, position, energy and spin

    write(*,('(25x,a14,6x,a9)')) 'Time iteration','Real time'
    do t_index=2,tmax                      !General loop over Time
       write(*,('(a20,8x,i4,10x,f11.6,x,a2)')) 'It is running-->Time',t_index,(t_index-1)*dt,'fs'
       starttime = omp_get_wtime()
       timem1=t_index-1                        !Time variables

       tscale=funcmod(t_index)
       tscale1=funcmod(timem1)

       !$OMP PARALLEL 
       !$OMP DO 
       do zposprim=1,zmax
          do E_index=1,Emax
             do S_index=1,2
                call Seffevaluat(S_index,E_index,zposprim,timem1) !This if is necessary to not repeat the Seffevaluat routine
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

       !do zposprim=1,zmax
       !   do S_index=1,2                
       !      call Seffvacuum_rut(zposprim,timem1,S_index)
       !   end do
       !end do

       do z_index=1,zmax                   !General loop over Positions
          iint%zpos=z_index 
          iint%time=t_index
          do S_index=1,2                   !General loop over Spin
             iint%sigma=S_index
             do E_index=1,Emax             !General loop over Energies
                iint%energy=E_index
                if(t_index==2)then
                   n(S_index,E_index,z_index,2)=Sexternal(S_index,E_index,z_index,2) 
                   WRITE (20+(S_index-1)*20+E_index,'(F28.8,$)') n(iint%sigma,iint%energy,iint%zpos,t_index)
                   flux(S_index,E_index,z_index,1)=0.0
                   if(verbosity == 2) call checkitems(0)
                   !cycle
                end if
                Rtau=tau(S_index,E_index,z_index,timem1)
                Item=exp(-dt/Rtau)*n(S_index,E_index,z_index,tscale1)+Seff(S_index,E_index,z_index,tscale1)
                z_ref=>iint%zpos    !I use a pointer here
                !
                value_phi=0.0               
                delta_t=0.0
                delta_t_tau=0.0
                sc_right=0.0
                sc_left=0.0
                !p_left=>reflect%z_sn(z_ref) !Pointer to  reflect%z_sn(z_ref) 
                flux_value=0.0
                !_________________________________________________________________________________________________
                !In the ruotine phi_right and phi_left "z_ref" represents the interfaces, NOT the spatial position
                !Therefore, they represent the amount of electrons crossing the interface "z_ref" from the right 
                ! and the left, respectively.
                !_________________________________________________________________________________________________
                if(z_ref.ne.1)then
                   call phi_right(1,z_ref,timem1,delta_t,delta_t_tau,value_phi)
                   Item=Item+value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref))
                   flux_value=value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref))                  !To compute the flux at any interface
                   spin_current_right(S_index,z_index,t_index)=spin_current_right(S_index,z_index,t_index)+value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref))/2.0
                   sc_right=sc_right+value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref))/2.0
                   value_phi=0.0
                end if
                call phi_left(2,z_ref,timem1,delta_t,delta_t_tau,value_phi)
                Item=Item+value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref))
                flux_value=flux_value+value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref))           !To compute the flux at any interface
                spin_current_left(S_index,z_index,t_index)=spin_current_left(S_index,z_index,t_index)+value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref))/2.0
                sc_left=sc_left+value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref))/2.0
                value_phi=0.0

                call phi_right(1,z_ref+1,timem1,delta_t,delta_t_tau,value_phi)
                Item=Item-value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref+1))
                flux_value=flux_value-value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref+1))         !To compute the flux at any interface
                spin_current_right(S_index,z_index,t_index)=spin_current_right(S_index,z_index,t_index)+value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref+1))/2.0
                sc_right=sc_right+value_phi*(1.0-reflect%ref_right(iint%sigma,iint%energy,z_ref+1))/2.0
                value_phi=0.0
                if(z_ref.ne.zmax)then
                   call phi_left(2,z_ref+1,timem1,delta_t,delta_t_tau,value_phi)
                   Item=Item-value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref+1))
                   flux_value=flux_value-value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref+1))      !To compute the flux at any interface
                   spin_current_left(S_index,z_index,t_index)=spin_current_left(S_index,z_index,t_index)+value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref+1))/2.0
                   sc_left=sc_left+value_phi*(1.0-reflect%ref_left(iint%sigma,iint%energy,z_ref+1))/2.0
                end if

                n(S_index,E_index,z_index,tscale)=Item
                WRITE (20+(S_index-1)*20+E_index,'(F28.8,$)') Item
                flux(S_index,E_index,z_index,tscale1)=flux_value 
                if(verbosity == 2) call checkitems(0)
                call SC_energy_dependent(.FALSE.,sc_right,sc_left)
             end do                          !Close general loop over Energies
             !Compute the number of electron in the Marco's sea
             MS=0.0
             do E_index=1,Emax
                MS=MS-Sexternal(S_index,E_index,z_index,t_index)
                n_total(S_index,z_index,t_index)=n_total(S_index,z_index,t_index)+n(S_index,E_index,z_index,tscale)
                flux_total(S_index,z_index,t_index)=flux_total(S_index,z_index,t_index)+flux(S_index,E_index,z_index,tscale1)
             end do
             nvacuum(S_index,z_index,t_index)=nvacuum(S_index,z_index,t_index-1)+Seffvacuum(z_index,t_index-1,S_index)+MS  
             write(97+S_index,'(F28.8,$)') nvacuum(S_index,z_index,t_index)
             write(620+S_index,'(F18.8,$)') n_total(S_index,z_index,t_index)
             write(622+S_index,'(F18.8,$)') flux_total(S_index,z_index,t_index)
             write(720+S_index,'(F18.8,$)') spin_current_right(S_index,z_index,t_index)
             write(722+S_index,'(F18.8,$)') spin_current_left(S_index,z_index,t_index)
          end do                        !Close general loop over Spin
          magnetization(z_index,t_index)=n_total(2,z_index,t_index)+nvacuum(2,z_index,t_index)-n_total(1,z_index,t_index)-nvacuum(1,z_index,t_index)
          flux_F(z_index,t_index)=flux_total(2,z_index,t_index)-flux_total(1,z_index,t_index)
          write(625,'(F18.8,$)') magnetization(z_index,t_index)
          write(626,'(F18.8,$)') flux_F(z_index,t_index)
       end do                     !Close general loop over Positions
       write(625,*) 
       write(626,*) 

       if(verbosity == 2) call checkitems(1)
       do S_index=1,2
          do E_index=1,Emax
             WRITE (20+(S_index-1)*20+E_index,*) 
          end do
          write(97+S_index,*) 
          write(620+S_index,*) 
          write(622+S_index,*) 
          write(720+S_index,*) 
          write(722+S_index,*) 
       end do
       endtime = omp_get_wtime()
       write(*,*), endtime-starttime
    end do                              !Close general loop over Time
    deallocate(nvacuum)


    write(*,*) '========================================================================='
    write(*,'(10x,a60)') 'Output files print in created directory: '//directory//''
    write(*,'(30x,l2)') succesful
    write(*,*) '========================================================================='
    write(*,'(25x,a10)') 'Total Time'
    write(*,'(16x,a3,25x,a3)') 'Sec','Min'
    write(*,'(6x,f20.15,8x,f20.15)') endtime-verystarttime,(endtime-verystarttime)/60.d0
    write(*,*) '========================================================================='
    write(*,'(a30,3x,a35)') 'Number of threads used','Number of processor in the computer'
    write(*,'(15x,i3,25x,i3)') OMP_get_max_threads(),OMP_get_num_procs()
    write(*,*) '========================================================================='
  end program interface








    





