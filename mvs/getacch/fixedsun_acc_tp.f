c---------------------------------------------------------------------
c*************************************************************************
c                           FIXEDSUN_ACC_TP.F
c*************************************************************************
c FIXEDSUN_ACC_TP returns the x,y,z components of the acceleration
c on NTP test particles due to the Sun on a fixed orbit about the system
c*************************************************************************
c  Input:
c
c  ntp                  ==>  number of test particles (integer)
c
c  istat                ==>  status of the test paricles (integer array)
c
c  xht(*),yht(*),zht(*) ==>  heliocentric positions of test particles 
c                             (real(8) vectors)
c************************************************************************* 
c      Output:
c
c  a_sun_x(*), a_sun_y(*), a_sun_z(*)  ==>  heliocentric components of 
c                                               acceleration (real(8) vectors) 
c*************************************************************************
c      Common:
c  msun                 ==>  mass of sun
c
c  ahel                 ==>  semimajor axis of orbit (assumes 0 i/e)
c
c*************************************************************************
c      Remarks:  
c      Author :  Nathan Kaib
c      Date   :  02/07/18
c      Last revision: 02/07/18 by Ramon Brasser
c*************************************************************************
      subroutine fixedsun_acc_tp(ntp,istat,xht, yht, zht,a_sun_x,
     &     a_sun_y,a_sun_z,time,mass)

       include '../../swift.inc'

c...  Inputs:
      integer ntp, istat(NTPMAX,NSTAT)
c
      real*8 xht(NTPMAX), yht(NTPMAX), zht(NTPMAX)
      real*8 mass(NPLMAX)

c...  Output
      real*8 a_sun_x(NTPMAX), a_sun_y(NTPMAX), a_sun_z(NTPMAX)

c...  Common
      real*8 time
      real*8 msun,ahel

c...  Internals
      integer i
      real*8 mcent
      real*8 perhel,angorb
      real*8 xsun,ysun,zsun
      real*8 rsun,rsun2,rsun3
      real*8 dist,dist2,dist3
      real*8 dx,dy,dz
      real*8 aprimx,aprimy,aprimz

      common/c_sunacc/msun,ahel
c---
c...  executable code
      mcent = mass(1) + msun
      perhel = 2.d0 * PI * sqrt(ahel**3.d0  / mcent)
      angorb = time/perhel * 2.d0 * PI
      xsun = ahel * cos(angorb)
      ysun = ahel * sin(angorb)
      zsun = 0.d0

      rsun2 = xsun*xsun + ysun*ysun + zsun*zsun
      rsun = sqrt(rsun2)
      rsun3 = rsun2 * rsun

      aprimx = mcent * xsun / rsun3
      aprimy = mcent * ysun / rsun3
      aprimz = mcent * zsun / rsun3
      do i = 1, ntp
         if(istat(i,1) .eq. 0) then
	   dx = xht(i) - xsun
	   dy = yht(i) - ysun
	   dz = zht(i) - zsun

	   dist2 = dx*dx + dy*dy + dz*dz
	   dist = sqrt(dist2)
	   dist3 = dist2 * dist

           a_sun_x(i) = - mcent * dx / dist3 - aprimx
           a_sun_y(i) = - mcent * dy / dist3 - aprimy
           a_sun_z(i) = - mcent * dz / dist3 - aprimz
         endif
      enddo

      return         
      end                       !  fixedsun_acc_tp.f

