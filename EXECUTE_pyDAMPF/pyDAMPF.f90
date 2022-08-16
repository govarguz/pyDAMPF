
!=======================pyDAMPF==============================
!The code written below is part of the pyDAMPF project and has
!been developed as a joint project between the IJS and UMSA. This 
!software will enable students and researchers from the 
!AFM community to easily plan AFM experiments.
!======================CREDITS==============================
!
! PI Horacio V. Guzman Ph.D.
! B.Sc. Willy Menacho N.
!
!We thank you in advance for sending your feedback and/or 
!suggestions to:
!             horacio.v.g@gmail.com
!
!======================pyDAMPF LICENSE==============================
! 
!Copyright (C) 2022  Horacio V. Guzman and Willy Menacho N.
! 
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
! 
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
! 
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
! 
! 
! What we observe is not nature itself, but nature exposed to
! our method of questioning. 
!                                  W. Heisenberg
!###########################################################
!====================================================================


!====================POINT MASS MODEL=======================
MODULE mypmmoda
  IMPLICIT NONE
  REAL*8, PARAMETER :: pi=3.141592654
  REAL*8, PARAMETER :: kb=1.380658e-23
  REAL*8, PARAMETER :: mu0=4*pi*1.e-7
  INTEGER, PARAMETER :: nmax=10 !diff
  INTEGER, PARAMETER :: ndin=12800
  REAL*8,DIMENSION(nmax) :: dydx,y,yout
  REAL*8 :: xminimo,xmaximo,dminimo,vminimo,intermax,amplitud,zcmin,dzc,zcmax,interzc,factor
  REAL*8 :: saltot,omega,omega0,periodo,kc,A0,F0,Q,rad,fase1,amp1,defl,virial1,pts1,ets1,fmedia,radSam
  REAL*8 :: tmin,tmax,fthermal0,Temp
  REAL*8 :: t,zc,a00,inter
  REAL*8 :: ham,ebarra,epsilon1,eta,epsilon2,ljmin,length,mtip,msample,b00
  REAL*8 :: eps,eps0, sigmat,sigmas,debye             !WM
  REAL*8 :: vdw,Hertz,visco,cap,lj,velo,dlvo          !WM
  REAL*8 :: npp,nper,nperfin,aux 
  REAL*8 :: fixzc,tcs,m,d,dmaximo,inte	!new
  INTEGER :: n,nzc,j,l,i,jmin,idum,label,tcsum,flagin,nzcfix
  CHARACTER*30  nombre
  !#####
  !FUNCTIONS AND SUBROUTINES
  !#####
  CONTAINS
  !#####
  SUBROUTINE mainbim()
  IMPLICIT NONE
  INTEGER :: j,l,i,jj
  REAL*8,DIMENSION(ndin) :: tarray,xarray,varray,xaux,dtarray,ftarray
  REAL*8,DIMENSION(ndin) :: vdwtarray,Hertztarray,viscotarray,captarray,ljtarray,dlvotarray !wm
  PRINT *, ndin
  !------------------------------Using an existing subroutine--------------------------------
  CALL input
  !---------------------------------Using a new subroutine-----------------------------------
  !INITIAL CONDITIONS
  ftarray(1)=0.0	!new
  dtarray(1)=0.0	!new
  !wm
  vdwtarray(1)=0.0
  Hertztarray(1)=0.0
  viscotarray(1)=0.0
  captarray(1)=0.0
  ljtarray(1)=0.0
  dlvotarray(1)=0.0
  !wm
  !----------------new
  y(1)=0.	
  y(2)=0.
  xarray(1)=0.
  varray(1)=0.
  tarray(n)=tmin+(n-1)*saltot
  ! # # # #
  !=============================START ZC ANALYSIS=======================================
  OPEN(1,FILE='zcdom.dfo')
  OPEN(2,FILE='tdom.dfo')
  DO l=1,Nzc  !Empieza el bucle en zc
  zc=zcmax-(l-1)*dzc 
  tcsum=0
  !=============================START TIME ANALYSIS=======================================
  DO j=1,ndin-1	 
  t=tmin+(j-1)*saltot	
  CALL derivs(t,y,dydx,zc,inter,d)	!new
  ftarray(j+1)=inter	!new
  dtarray(j+1)=d		!new
  IF (dtarray(j+1)>0.0) THEN
  tcsum=tcsum+1		!new
  END IF
  CALL rk4(y,dydx,8,t,saltot,yout,zc)
  xarray(j+1)=yout(1)
  varray(j+1)=yout(2)
  tarray(j)=t
  !wm
  call viscoel(dtarray(j+1),varray(j+1),visco)
  viscotarray(j+1)=visco
  call capilar(dtarray(j+1),varray(j+1),cap)
  captarray(j+1)=cap
  call attractive(dtarray(j+1),varray(j+1),vdw)
  vdwtarray(j+1)=vdw
  call repulsive(dtarray(j+1),Hertz)
  Hertztarray(j+1)=Hertz
  call lej(dtarray(j+1),lj)
  ljtarray(j+1)=lj
  call dlvoe(dtarray(j+1),dlvo)
  dlvotarray(j+1)=dlvo
  !wm
  DO i=1,2
  y(i)=yout(i)
  END DO
  !=============================ENDS TIME ANALYSIS=======================================
  END DO
  nzcfix=int(1.99999+(fixzc-zcmin)/dzc)
  nzcfix=int(nzc-nzcfix)
  !=============================Start APD analysis===================================
  xmaximo=0.
  xminimo=0.
  virial1=0.
  pts1=0.
  fmedia=0.
  DO j=(nper-nperfin)*npp+1,ndin
  ! Máximo y mínimo		
  IF (xarray(j).lt.xminimo) THEN
  xminimo=xarray(j)
  jmin=j
  END IF
  IF (xarray(j).gt.xmaximo) THEN
  xmaximo=xarray(j)
  END IF		
  ! Virial y energía disipada
  CALL interaction(xarray(j)+zc,varray(j),inte)
  !extended simpson's rule (numerical recipes)										
  IF((j==(nper-nperfin)*npp).OR.(j==n)) THEN
  factor=1./3.
  ELSE IF ((j/2.)==int(j/2.)) THEN
  factor=4./3.
  ELSE IF ((j/2.)/=int(j/2.)) THEN
  factor=2./3.	
  END IF
  virial1=virial1+factor*xarray(j)*inte/(nperfin*npp)
  pts1=pts1+factor*varray(j)*inte/(nperfin*npp)
  fmedia=fmedia+factor*inte/(nperfin*npp)
  END DO
  amplitud=(xmaximo-xminimo)/2.
  defl=(xmaximo+xminimo)/2.	
  vminimo=varray(jmin)
  dminimo=zc+xminimo
  dmaximo=zc+xmaximo
  tcs=(1-(tcsum*1.0/(ndin-1)*1.0))
  CALL interaction(dminimo,vminimo,intermax)
  ets1=pts1*periodo														   		
  ! Transformada de Fourier.
  amp1=0.
  fase1=0.
  ! call ampphase(xarray,tarray,omega,ndin,fase1,amp1)

      !Escribir al fichero
  WRITE(1,1000) zc*1.e9,amp1*1.e9,fase1,dminimo*1.e9,dmaximo*1.e9,intermax*1.e9,-ets1*1.e20&
  &,virial1*1.e20,fmedia*1.e9,defl*1.e9,tcs,pts1*1.e16
  !Maybe the if sentence for time domain must come here instead of going where it was before.
  nzcfix=int(1.99999+(fixzc-zcmin)/dzc)
  nzcfix=int(nzc-nzcfix)
  !=============================T ANALYSIS PRINTING=======================================
  !write(1,2001) 't/T','zt','vt','z1','v1','z2','z2','z3','v3','z4','v4','zc','Finst'
  !write(1,2001) 'adim','nm','m/s','nm','m/s','nm','m/s','nm','m/s','nm','m/s','nm 
  !write (1,2000) tarray(j)/periodo, xarray(j)*1.e9,varray(j),zc*1.e9
  !changed tolerance
  IF (abs(l-nzcfix).LE.0.1) THEN
  DO jj=(nper-nperfin)*npp+1,ndin
  ! PRINT*,'This condition is working'
  WRITE(2,2000) tarray(jj)/periodo,xarray(jj)*1.e9,varray(jj),zc*1.e9,ftarray(jj)*1.e9&
  &,dtarray(jj)*1.e9,vdwtarray(jj)*1.e9,Hertztarray(jj)*1.e9,viscotarray(jj)*1.e9&
  &,captarray(jj)*1.e9,ljtarray(jj)*1.e9,dlvotarray(jj)*1.e9
  END DO
  END IF
  END DO
  CLOSE(1)
  CLOSE(2)
  1000 FORMAT(12F20.5)
  2000 FORMAT(12F20.5)
  END SUBROUTINE
  !# # #
  !END OF MAIN SUBROUTINE
  !# # #
  SUBROUTINE rk4(y,dydx,n,x,h,yout,zc)
  IMPLICIT NONE
  INTEGER :: i,n,nmax
  PARAMETER (nmax=50)
  REAL*8 h,x,dydx(nmax),y(nmax),yout(nmax)
  REAL*8 h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax),zc
  !BEGINS
  hh=h*0.5
  h6=h/6.
  xh=x+hh
  DO 11 i=1,n
  yt(i)=y(i)+hh*dydx(i)
  11    continue
  call derivs(xh,yt,dyt,zc,inter,d)
  DO 12 i=1,n
  yt(i)=y(i)+hh*dyt(i)
  12    continue
  call derivs(xh,yt,dym,zc,inter,d)
  DO 13 i=1,n
  yt(i)=y(i)+h*dym(i)
  dym(i)=dyt(i)+dym(i)
  13    continue
  call derivs(x+h,yt,dyt,zc,inter,d)
  DO 14 i=1,n
  yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
  14    continue
  RETURN
  END SUBROUTINE
  !#####  
  SUBROUTINE derivs(x,y,dydx,zc,inter,d)
  IMPLICIT NONE
  INTEGER :: i,nmax
  PARAMETER (nmax=50)
  REAL*8 :: dydx(nmax),y(nmax)
  REAL*8 :: zc,d,inter,veloci,x
  INTENT(IN) x,y,zc
  INTENT(OUT) dydx,d,inter
  !BEGINS
  d=zc+y(1)
  veloci=y(2)
  CALL interaction(d,veloci,inter)
  dydx(1)=y(2)
  dydx(2)=-y(1)*omega0**2-omega0*y(2)/Q+(F0*DCOS(omega0*x)+inter)/m
  END SUBROUTINE		

  !#####
  SUBROUTINE viscoel(d,velo,visco)
  IMPLICIT NONE
  INTENT(IN) d,velo
  INTENT(OUT) visco
  REAL*8 :: d,velo,visco
  IF(d<a00) THEN
  visco=-eta*dsqrt(rad*(a00-d))*velo		!¿No sería a00-d en lugar de dabs(d)?
  ELSE IF(d>=a00) THEN
  visco=0.
  END IF
  END SUBROUTINE
  !#####
  SUBROUTINE capilar(d,velo,cap)
  IMPLICIT NONE
  INTENT(IN) d,velo
  INTENT(OUT) cap
  REAL*8 :: d,velo
  REAL*8 :: cap
  IF (d>13.8*a00) THEN
  cap=0.
  ELSE IF ((d>a00).AND.(d<13.8*a00)) THEN 
  cap=-0.5*epsilon2*Rad*(1.+velo/abs(velo))*(1.-d/(13.8*a00))		!Dc=2.27 nm=13.8*a00
  ELSE IF	(d<a00) THEN
  cap=0.
  END IF
  END SUBROUTINE
  !#####
  SUBROUTINE repulsive(d,Hertz)
  IMPLICIT NONE
  INTENT(IN) d
  INTENT(OUT) Hertz
  REAL*8 :: d
  REAL*8 :: Hertz 
  if ((d)>a00) then
  Hertz=0.
  else if (d<a00) then
  Hertz=(Ebarra*sqrt(rad)*(a00-d)**1.5)
  !Hertz=0.
  end if
  END SUBROUTINE
  !#####
  SUBROUTINE attractive(d,velo,vdw)
  IMPLICIT NONE
  INTENT(in) d,velo
  INTENT(out) vdw
  REAL*8 :: d, velo
  REAL*8 :: vdw
  if (d>a00) then
  vdw=-Ham*Rad/(6.*d**2)						!-(Ham*Rad/(6.*d**2))*(1+(d/Rad)/(2+d/Rad))**2
  else if (d<a00) then
  vdw=-Ham*Rad/(6.*a00**2)
  end if
  if(abs(velo).gt.1.e-20) then
  vdw=vdw*(1.+epsilon1*velo/abs(velo))
  end if	

  END SUBROUTINE
  !########DLVO###########    !WM
  SUBROUTINE dlvoe(d,dlvo)
  IMPLICIT NONE
  INTENT(in) d
  INTENT(out) dlvo
  REAL*8 :: d
  REAL*8 :: dlvo
  if (d>a00) then
  dlvo=(((4.*pi*Rad)/(eps*eps0))*sigmat*sigmas*debye*(EXP((-d/debye))))
  else if (d<a00) then
  dlvo=(((4.*pi*Rad)/(eps*eps0))*sigmat*sigmas*debye*(EXP((-a00/debye))))
  end if
  END SUBROUTINE
  !##### lennard-jones
  SUBROUTINE lej(d,lj)
  IMPLICIT NONE
  INTENT(IN) d
  INTENT(OUT) lj
  REAL*8 :: d
  REAL*8 :: lj,ljaux
  IF (d>a00) THEN
  lj=1.5*ljmin*(length**2)*(length**4/(3.*d**6)-1./d**2)      !dForce old
  ELSE IF (d<a00) THEN
  lj=0.
  END IF 
  END SUBROUTINE            !WM
  !##### magnetic
  SUBROUTINE magnetic(d,mag)
  IMPLICIT NONE
  INTENT(IN) d
  INTENT(OUT) mag
  REAL*8 :: d
  REAL*8 :: mag,delta
  !ferritina antigua
  !mag_width=10.e-9	!=radio de las nanoparticulas (= d_m en Rasa et al. JColl Int Sci 2002)
  !mag_thick=0.8e-9	!=grosor de la capa magnetica en la punta (= delta en Rasa et al. JColl Int Sci 2002)
  !mag=-mu0*(8.*pi**2/3.)*Mtip*Msample*(mag_width*(rad-mag_thick))**3*(d+0.5*rad+0.5*mag_width+mag_thick)**(-4)
  !ferritina nueva
  !Mtip y Msample son momentos magnéticos, y no magnetizaciones (=momento magnético por unidad de volumen).
  delta=12.e-9 ! es el radio total de la ferritina
  mag=-(1.5*mu0/pi)*Mtip*Msample/(rad+delta+d)**4
  END SUBROUTINE	
  !##### 
  SUBROUTINE interaction(d,velo,inter)
  IMPLICIT NONE
  INTENT(IN) d,velo
  INTENT(OUT) inter
  REAL*8 :: d,velo
  REAL*8 :: inter,vdw,Hertz,visco,cap,lj,mag,dlvo
  call repulsive(d,Hertz)
  call attractive(d,velo,vdw)
  call viscoel(d,velo,visco)
  call capilar(d,velo,cap)
  call dlvoe(d,dlvo)
  call lej(d,lj)
  call magnetic(d,mag)
  inter=vdw+Hertz+visco+cap+lj+mag+dlvo
  END SUBROUTINE		  
  !##### 
  ! SUBROUTINE ampphase(xarray,tarray,omega,ndin,fase,amp)
  ! IMPLICIT NONE
  ! INTENT(IN) :: xarray,tarray,omega,ndin
  ! INTENT(OUT) :: fase,amp
  ! REAL*8 :: fase,amp,omega
  ! COMPLEX(8) fourier
  ! COMPLEX(8),PARAMETER:: i=(0.,1.)
  ! REAL*8 :: xarray(ndin),tarray(ndin)
  ! INTEGER :: ndin,j
  ! fourier=0.
  ! do j=(nper-nperfin)*npp+1,ndin
  ! fourier=fourier+xarray(j)*cdexp(-i*omega*tarray(j))  
  ! end do
  ! amp=2.*cdabs(fourier)/(nperfin*npp)
  ! fase=DACOS(dreal(fourier)/cdabs(fourier))*(180/pi)
  ! END SUBROUTINE
  !##### 	
  SUBROUTINE input
  IMPLICIT NONE
  CHARACTER(16) colum
  OPEN(10,file='inputs.ini')
  READ(10,*) colum,omega0
  READ(10,*) colum,omega
  READ(10,*) colum,kc
  READ(10,*) colum,Q
  !READ(10,*) colum,Temp
  READ(10,*) colum,A0
  READ(10,*) colum,rad
  READ(10,*) colum,a00
  READ(10,*) colum,Ebarra
  READ(10,*) colum,Ham
  READ(10,*) colum,eta
  READ(10,*) colum,epsilon1
  READ(10,*) colum,epsilon2
  READ(10,*) colum,ljmin
  READ(10,*) colum,length
  READ(10,*) colum,Mtip
  READ(10,*) colum,Msample
  READ(10,*) colum,Nper
  READ(10,*) colum,Npp
  READ(10,*) colum,Nperfin
  READ(10,*) colum,zcmax
  READ(10,*) colum,zcmin
  READ(10,*) colum,dzc
  READ(10,*) colum,fixzc
  READ(10,*) colum,eps
  READ(10,*) colum,eps0
  READ(10,*) colum,sigmat
  READ(10,*) colum,sigmas
  READ(10,*) colum,debye
  CLOSE(10)
  idum=-1
  n=nper*npp
  F0=kc*A0/Q
  omega0=omega0*2.*pi
  omega=omega*2.*pi
  periodo=2.*pi/omega
  tmax=nper*periodo;
  tmin=0.
  saltot=(tmax-tmin)/n
  nzc=int(1.99999+(zcmax-zcmin)/dzc)
  m=kc/omega0**2
  END SUBROUTINE
	
END MODULE mypmmoda
