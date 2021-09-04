  program adsorption

    implicit none
    real*8 :: Izero,Pw_total,z_pos_i,z_pos_f
    real*8,allocatable,dimension(:) :: refrac_index,thickness,alpha,volume,power,intensity,extinc_coeff,epsilon_2
    integer :: num_lay,i

    !Open and read input variables
    open (UNIT = 1, FILE = "input_adsorption2", STATUS = "OLD", ACTION = "READ")

    read(1,*) !Number of different layers
    read(1,*) num_lay

    allocate(refrac_index(num_lay),extinc_coeff(num_lay),thickness(num_lay),alpha(num_lay),volume(num_lay))
    allocate(power(num_lay),intensity(num_lay),epsilon_2(num_lay))

    read(1,*) !refrac_index
    read(1,*) refrac_index(1:num_lay)
    read(1,*) !Extinction coefficient
    read(1,*) extinc_coeff(1:num_lay)
    read(1,*) !thickness  (nm)
    read(1,*) thickness(1:num_lay)
    read(1,*) !alpha (nm^{-1})
    read(1,*) alpha(1:num_lay)
    read(1,*) !Volume (cm^3/mol)
    read(1,*) volume(1:num_lay)
    read(1,*) !Izero (10^15 eV/cm^2)
    read(1,*) Izero

    !Calculate the adsorbed power
    power(1)=refrac_index(1)*(exp(-thickness(1)*alpha(1))-1.d0)
    Pw_total=power(1)
    z_pos_i=thickness(1)
    do i=2,num_lay
       z_pos_f=z_pos_i+thickness(i)
       write(*,*) z_pos_i,z_pos_f
       power(i)=refrac_index(i)*(exp(-z_pos_f*alpha(i))-exp(-z_pos_i*alpha(i)))
       Pw_total=Pw_total+power(i)
       z_pos_i=z_pos_f
    end do
    do i=1,num_lay
       intensity(i)=Izero*power(i)*volume(i)/(Pw_total*thickness(i)*6.022*97.5d0)
       write(*,('(3x,f6.3,x,a1,3x,f9.5)')) power(i)/Pw_total*100.d0,"%",intensity(i)
    end do
    
    !Calculation for homogeneous adsorption
    write(*,*) "   Homogeneous adsorption   "
    Pw_total=0.0
    do i=1,num_lay
       epsilon_2(i)=2.0*refrac_index(i)*extinc_coeff(i)
       Pw_total=Pw_total+thickness(i)*epsilon_2(i)
    end do
    do i=1,num_lay
       power(i)=thickness(i)*epsilon_2(i)/Pw_total
       intensity(i)=Izero*power(i)*volume(i)/(thickness(i)*6.022*97.5d0)
       write(*,('(3x,f6.3,x,a1,3x,f9.5)')) power(i)*100.d0,"%",intensity(i)
    end do

  end program adsorption



