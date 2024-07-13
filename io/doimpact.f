       subroutine doimpact(ntp,istat,mass,vxht,vyht,vzht,mimp,ximp,
     &           yimp,zimp,vximp,vyimp,vzimp,
     &           impobj)

       include '../swift.inc'
       include '../stars2.inc'

       integer ntp
       integer istat(NTPMAX,NSTAT)
       real*8 mass(NPLMAX)
       real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
       real*8 mimp,ximp,yimp,zimp,vximp,vyimp,vzimp
       integer impobj

       real*8 rdotv,costhet,rmag,rmag2,vmag,vmag2
       real*8 delvx,delvy,delvz,massrat

c-----
c...  Executable code 

       if (istat(impobj,1).eq.0) then
	  rmag2 = ximp * ximp +  yimp * yimp +  zimp * zimp
	  rmag = sqrt(rmag2)

	  vmag2 = vximp * vximp +  vyimp * vyimp +  vzimp * vzimp
	  vmag = sqrt(vmag2)

          rdotv = ximp * vximp + yimp * vyimp + zimp * vzimp
	  costhet = rdotv / rmag / vmag

	  massrat = mimp / (mass(1) + mimp)
	  delvx = -massrat * costhet * vximp !see Eq 2 of Nesvorny 2019
	  delvy = -massrat * costhet * vyimp 
	  delvz = -massrat * costhet * vzimp 

	  vxht(impobj) = vxht(impobj) + delvx
	  vyht(impobj) = vyht(impobj) + delvy
	  vzht(impobj) = vzht(impobj) + delvz
       endif

       return
       end
