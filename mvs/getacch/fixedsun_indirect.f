c---------------------------------------------------------------------
c*************************************************************************
c                           FIXEDSUN_INDIRECT.F
c*************************************************************************
c FIXEDSUN_INDIRECT returns the x,y,z components of the indirect acceleration
c simulation particles due to the Sun on a fixed orbit about the system
c*************************************************************************
c  Input:
c
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
c      Date   :  01/28/19
c*************************************************************************
      subroutine fixedsun_indirect(a_sun_x,
     &     a_sun_y,a_sun_z,time,mass)

       include '../../swift.inc'

c...  Inputs:
      integer ntp, istat(NTPMAX,NSTAT)
c
      real*8 mass(NPLMAX)

c...  Output
      real*8 a_sun_x, a_sun_y, a_sun_z

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

      a_sun_x = 0.0d0 - aprimx
      a_sun_y = 0.0d0 - aprimy
      a_sun_z = 0.0d0 - aprimz

      return         
      end                       !  fixedsun_acc_tp.f

