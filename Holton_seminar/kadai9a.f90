program kadai9a

 implicit none

 !変数の指定
 integer ,parameter :: cp=1000 , s=1367
 real ,parameter ::  rho=1.2 , h=8300 ,  e=0.6 , sig=0.0000000567, dm=3.0 !month interval (change)
 real ,parameter :: dt=2592000 ! time interval change(second)(1month)
 integer,parameter :: nmax = 1000 !time step (nochange)
 real :: a,k1,k2,k3,k4
 real :: T(0:nmax-1)
 integer :: n

 !出力ファイル
 open(49,file='kadai9a_1-5.dat',status='replace',form='formatted')
 open(50,file='kadai9a_2-5.dat',status='replace',form='formatted')

 T(0)=293.0
 !計算(euler)

 do n= 0,nmax-1

     !アルベドの条件わけ
     if ( T(n) - 273.15 > -10 )then
     a=0.30
     else
     a=0.62
     endif

    T(n+1) = T(n) + ((s*(1-a))/(4*cp*rho*h)-(e*sig*(T(n))**4)/(cp*rho*h))*dt*dm


   write(49,'(f5.2,2x,f15.7)') n*dm,T(n)

  ! if( abs( T(n+1) - T(n) ) < 0.001 ) exit

 end do

 close(49)

 !計算(lunge-kutta)

 do n= 0,nmax-1

   !アルベドの条件わけ
  if ( T(n) - 273.15 > -10)then
  a=0.30
  else
  a=0.62
  endif

   k1=(s*(1-a))/(4*cp*rho*h)-(e*sig*(T(n))**4)/(cp*rho*h)
   k2=(s*(1-a))/(4*cp*rho*h)-(e*sig*(T(n)+(dt*dm*k1/2))**4)/(cp*rho*h)
   k3=(s*(1-a))/(4*cp*rho*h)-(e*sig*(T(n)+(dt*dm*k2/2))**4)/(cp*rho*h)
   k4=(s*(1-a))/(4*cp*rho*h)-(e*sig*(T(n)+(dt*dm*k3))**4)/(cp*rho*h)

   T(n+1) = T(n) + (dt*dm/6)*(k1+2*(k2+k3)+k4)
   write(*,*)T(n)

 write(50,'(f5.2,2x,f12.5)') n*dm , T(n)

 !if( abs( T(n+1) - T(n) ) < 0.001 ) exit

 enddo

 close(50)

end program kadai9a
