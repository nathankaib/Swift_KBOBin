c*************************************************************************
c                            RMVS4_STEP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c both massive and test particles
c INCLUDES close approuches between test particles and planets
c
c VERSION 4 of RMVS.  The only difference from RMVS3 is that the forces are 
c                     always calculated between the planets.  This hopefully
c                     will fix the problem that when I do multiple runs, the 
c                     planets end up at different places.
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial planet position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial planet velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial tp position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final planet position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final planet velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c
c
c Remarks: Based on rmvs3_step.f
c Authors:  Hal Levison 
c Date:    9/18/05
c Last revision: 

      subroutine rmvs4_step(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt)	

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer i1sttp,i1stpl,i1sto,icflg,i,j,np
      real*8 xtmp(NPLMAX,NTENC),ytmp(NPLMAX,NTENC)
      real*8 ztmp(NPLMAX,NTENC),rts,peri(NTPMAX)
      real*8 vxtmp(NPLMAX,NTENC),vytmp(NPLMAX,NTENC)
      real*8 vztmp(NPLMAX,NTENC)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)
      Integer nenc(NPLMAX),itpenc(NTPMAX,NPLMAX),ienc(NTPMAX)
      integer istattmp(NTPMAX,NSTAT),isperi(NTPMAX)

c----
c...  Executable code 

      i1sttp = i1st
      i1sto = i1st

c...  Note: were not doing close encounters since we just have the Sun and Primary,
c...  and the Sun is on a fixed orbit
      call rmvs4_step_norm(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,xh,yh,
     &     zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,dt)

      do i=1,ntp
 	if(istat(i,1).eq.0) then
 	   istat(i,2) = 0
 	endif
      enddo

      return       !  NOTE AN EXIT
c      endif

      end   ! step_enc
c------------------------------------------------------------------------






