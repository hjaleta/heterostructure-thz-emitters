module functions
  use omp_lib

  implicit none
  !integer, parameter :: ntpnt=1000
  public B5, trapez

     CONTAINS
        
       !Rutina B5.F
       !================================================================
       !JB=1,2,3 para integrar en 1,2,3 bloques a partir del punto de
       !abcisas IAS.  
       !M1,M2,M3 numero de pasos en cada bloque, H1,H2,H3
       !paso correspondiente. HI=HI/22.5 
       !================================================================
       
       real FUNCTION B5(JB,M1,M2,M3,H1,H2,H3,A,IAS)
         IMPLICIT NONE
         !integer, parameter :: ntpnt=1000
         real :: HL5,H5,S2,S3,S4,AX
         integer :: mp,j0,ir,N,nb,n1,L,I,mp1
         integer,intent(in) :: JB,m1,m2,m3,ias
         real,intent(in) :: h1,h2,h3
         real,dimension(:) :: A
         
         B5=0.D0
         N=M1-IAS+1
         MP=M1
         H5=H1
         J0=IAS
         IR=JB-2
1        if(N.le.0) goto 60
         if(N.lt.4) N=4
         IF(N-4) 44,2,2
2        HL5=H5/32.D0
         NB=N/4
         N1=N-4*NB+1
         S2=0.D0
         S3=S2
         S4=S2
         L=J0
         DO 4 I=1,NB
            S2=S2+A(L+2)
            S3=S3+A(L+1)+A(L+3)
            S4=S4+A(L)+A(L+4)
4        L=L+4
            GOTO (8,12,16,20),N1
8        if(MP.lt.4) then
               L=MP+J0
               MP1=MP
         else
               MP1=4
         end if
         GOTO (20,16,12,10),mp1
10       AX=0.d0
         GOTO 40
12       AX=HL5*(-19.D0*A(L-3)+106.D0*A(L-2)-264.D0*A(L-1)+646.D0*A(L)+&
              251.D0*A(L+1))
         GOTO 40
16       AX=HL5*(-8.D0*A(L-2)+32.D0*A(L-1)+192.D0*A(L)+992.D0*A(L+1)+&
              232.D0*A(L+2))
         GOTO 40
20       AX=HL5*(-27.D0*A(L-1)+378.D0*A(L)+648.D0*A(L+1)+918.D0*A(L+2)+&
              243.D0*A(L+3))
40       if(mp.lt.4) then
            B5=B5+H5*(7.D0*S4+32.D0*S3+12.D0*S2)-AX
         else
            B5=B5+H5*(7.D0*S4+32.D0*S3+12.D0*S2)+AX
         end if
         IF(IR) 60,50,55
50       J0=M1+IAS
         N=M2
         MP=M2
         H5=H2
         IR=-1
         GOTO 1
55       J0=M2+M1+IAS
         N=M3
         H5=H3
         MP=M3
         IR=0
         GOTO 1
60       RETURN
          
44       WRITE(4,200)
200      FORMAT(//1X,'  NUMERO DE INTERVALOS MENOR QUE 4'/)
         STOP
       END FUNCTION B5     


       !Function to integrate the function A using the Trapezoidal Method      
       !         I(A)=hstep/2*(A(1)+2*sum_{i=1}^{minter-1} A(i)+A(minter))
       Real function trapez(hstep,minter,A)
         implicit none
         real, intent(in) :: hstep
         integer :: i
         integer, intent(in) :: minter
         real,dimension(:) :: A
         
         trapez=0.0
         do i=1,minter
            if(i.ne.1.and.i.ne.minter)then
               trapez=trapez+hstep*A(i)
            else
               trapez=trapez+hstep*A(i)/2.0
            end if
         end do
         
       end function trapez


       !Function to integrate the function A using the Trapezoidal Method      
       !         I(A)=hstep/2*(A(1)+2*sum_{i=1}^{minter-1} A(i)+A(minter))
       Real function rectangular(hstep,minter,A)
         implicit none
         real, intent(in) :: hstep
         integer :: i
         integer, intent(in) :: minter
         real, intent(in), dimension(:) :: A
         
         rectangular=0.0
         do i=1,minter
            rectangular=rectangular+hstep*A(i)
         end do
         
       end function rectangular
       
       !Function to evaluate the External source at a determinate space and time,
       !and for a given energy and spin. It depends on the material in which 
       !the signal is.
       real function Sexternal(sigma,energy,zpos,time)
         use global, only: nummedia, placeinterfaces, excavrel, &
              lambdacoefficients, timewidth, t0source, dt, dz, extsourcetype
         implicit none
         REAL,PARAMETER	::  pi=3.14159265
         integer :: zmin, zmax
         integer :: i,j
         integer, intent(in) :: sigma, energy, zpos, time
         !It should be convenient if we write the expresion for each case
         select case (extsourcetype)
         case('delta')       !if(extsourcetype=='delta') then
            !if(zpos.eq.10.and.time.eq.10) then
            if((zpos == 10.and.time == 10).or.(zpos == 20.and.time == 20).or.(zpos == 55.and.time == 152))then
               Sexternal=100.0
            !else if((zpos.eq.15.and.time.eq.11).or.(zpos == 81.and.time.eq.11)) then
            !   Sexternal=90.0
            else
               Sexternal=0.0
            end if
         case('test')       !else if(extsourcetype=='test') then
            Sexternal=10*exp(-zpos*dz/30.0-2.0*(time*dt-10.0)**2/(10.0**2))
         case('expgauss')       !else if(extsourcetype=='expgauss')then       
            do j=1,nummedia
               zmin=INT(placeinterfaces(j-1)/dz)+1
               zmax=INT(placeinterfaces(j)/dz)+1
               if(zpos.ge.zmin.and.zpos.lt.zmax)then
                  i=j
                  exit
               end if
            end do
!            Sexternal=(placeinterfaces(i)-placeinterfaces(i-1))*excavrel(sigma,energy,i)/&
!                 (lambdacoefficients(i)*timewidth*sqrt(pi*log(16.))*(exp(-placeinterfaces(i-1)/&
!                 lambdacoefficients(i))-exp(-placeinterfaces(i)/lambdacoefficients(i))))*&
!                 exp(-(zpos-0.5)*dz/lambdacoefficients(i))*exp(-4*log(2.0)*&
!                 ((time-1)*dt-t0source)**2/(timewidth**2))
!!            Sexternal=excavrel(sigma,energy,i)/(timewidth*sqrt(pi*log(16.)))*&
!!                 exp(-4*log(2.0)*((time-1)*dt-t0source)**2/(timewidth**2))
            Sexternal=excavrel(sigma,energy,i)*sqrt(log(16.))/(timewidth*sqrt(pi))*&
                 exp(-4*log(2.0)*((time-1)*dt-t0source)**2/(timewidth**2))
         case('expgaus2')       !else if(extsourcetype=='expgauss')then       
            do j=1,nummedia
               zmin=INT(placeinterfaces(j-1)/dz)+1
               zmax=INT(placeinterfaces(j)/dz)
               if(zpos.ge.zmin.and.zpos.le.zmax)then
                  i=j
                  exit
               end if
            end do
            Sexternal=excavrel(sigma,energy,i)*&
                 tanh(pi*timewidth/dt)*&
                 (exp(dz/lambdacoefficients(i))-1)/&
                 (exp(-(zmin-1)*dz/lambdacoefficients(i))- exp(- zmax*dz/lambdacoefficients(i)))
            Sexternal=Sexternal*exp(-zpos*dz/lambdacoefficients(i))* &
                 (dt/pi)*(timewidth/((time*dt-t0source)**2 + timewidth**2) )
         case('general')    !else if(extsourcetype=='general') then
            WRITE (*,*) 'Sexternal has been read throught out an input file'
         case default
            WRITE(*,*) 'Error in the select case'
         end select
         
         return
       end function Sexternal

       integer function funcmod(tprim)
         use global, only: tbound
         integer :: tprim
         
         funcmod=mod(tprim,tbound)
         if(funcmod.eq.0)funcmod=tbound

         return
       end function funcmod

       !Function to compute the electron probability of being scattered and to go to the Marco's sea
       !This is computed for every energy, spin and for every position.
       real function Seffvacuum(zpos,time,sigma)
         use global, only: Emax, n, tau, p, dE,dt
         implicit none
         real Ieprim, IEninteg
         real, dimension(:), allocatable :: QEninteger
         integer :: enerprim, sigprim, zpos,time, sigma, tscale
         
         allocate(QEninteger(Emax))
         tscale=funcmod(time)
         do enerprim=1,Emax
            QEninteger(enerprim)=0.0
            do sigprim=1,2
               QEninteger(enerprim)=QEninteger(enerprim)+p(sigprim,enerprim,sigma,Emax+1,zpos,time)&
                    *n(sigprim,enerprim,zpos,tscale)*(1.0-exp(-dt/tau(sigprim,enerprim,zpos,time)))
            end do
         end do
         !IEninteg=B5(1,Emax-1,0,0,dE/22.5,0.0,0.0,QEninteger,1)
         IEninteg=rectangular(dE,Emax,QEninteger)
         Seffvacuum=IEninteg
         deallocate(QEninteger)
         
         return
       end function Seffvacuum

  !Now we evaluate the next operator:
  !Seff(sigma,energy,zpos,time) Phi
  !where: Phi = eq(19) Battiato et al. J. App. Phys 115, 172611 (2014)
  !
  !
  !  life_timebot=Max(t-t0,life_time+dt)
  !  life_timetop=Max(t+dt-t0,life_time+dt)
  !
  !  taured=life_timetau+pscatt/(life_timetau+life_time)-1.0/tau
  !
  !We use the exponential integral Ei(x) to evaluate this quantity. It is calculated in the Logarithmic_Integral module.

       real function Phi(z_half,z0,t0,delta_t,delta_t_tau)
         use global, only: dt, tau, v, iint, Rtau,lftravel,pscatt
         use Logarithmic_Integral
         implicit none
         integer :: z,t,tprime
         integer, intent(in) :: z_half,z0,t0
         real :: t_m_deltat
         real, intent(in) :: delta_t,delta_t_tau
         real :: life_time, prob_scatt
         real :: R_t,R_t0
         real :: varmax, D_time_bot, D_time_top, tau_red
         real(dp) :: Qxtop, Qxbot

         t=iint%time-1                                 !Iteration time
         z=iint%zpos                                   !Iteration position
         R_t=(t-1)*dt                                  !Real value of integer t
         R_t0=(t0-1)*dt                                !Real value of integer t0
         t_m_deltat = R_t-abs(delta_t)

         tprime=int(t_m_deltat/dt)+1
         life_time=lftravel(z_half,z0,tprime,iint%zpos)                  !Calculation of life_time for + case. Real value   \equiv life_time(z+dz/2,t,z0)
         prob_scatt=pscatt(z_half,z0,tprime,iint%zpos)                   !Calculation of prob_scatt for + case. Real value

         varmax=R_t-R_t0
         D_time_bot=max(varmax,abs(delta_t+life_time)) !Real value for Delta_t_bot
         varmax=R_t+dt-R_t0
         D_time_top=max(varmax,abs(delta_t+life_time)) !Real value for Delta_t_top
         
         tau_red=(delta_t_tau+prob_scatt)/(delta_t+life_time)-1.0/Rtau ! Real value for tau_red
         if(D_time_bot == D_time_top)then
            Phi=0.0
         else
            Phi=(delta_t+life_time)/2.0*(exp(-(R_t+dt-R_t0)/Rtau-D_time_bot*tau_red)/D_time_bot- &
            exp(-(R_t+dt-R_t0)/Rtau-D_time_top*tau_red)/D_time_top)
            if(abs(tau_red).gt.0.001)then
               Qxtop=-D_time_top*tau_red
               Qxbot=-D_time_bot*tau_red
               Phi=Phi+(delta_t+life_time)/2.0*exp(-(R_t+dt-R_t0)/Rtau)*tau_red*(dei(Qxbot)-dei(Qxtop))
            end if
         end if
         return
       end function Phi

END MODULE functions

