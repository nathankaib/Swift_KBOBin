       subroutine kb_enc_add(nbod,mass,xh,yh,zh,vxh,vyh,vzh,rplsq,
     &           mkb,xkb,ykb,zkb,vxkb,vykb,vzkb,dt,impappflag,
     &           icenc,ntp,xht,yht,zht,vxht,vyht,vzht)

       include '../swift.inc'
       include '../stars2.inc'

       integer nbod,icenc,ntp
       real*8 mass(NPLMAX),dt
       real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
       real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
       real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
       real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
       real*8 r2hill(NPLMAX),vss,bmin,rplsq(NPLMAX)
       real*8 mkb,xkb,ykb,zkb,vxkb,vykb,vzkb

       real*8 rdotv,costhet,rmag

       integer impappflag

       common/c_rhill/r2hill

c-----
c...  Executable code 

       if (impappflag.eq.0) then
	   nbod = nbod + 1
	   icenc = icenc + 1

	   mass(nbod) = mkb
	   xh(nbod) = xkb
	   yh(nbod) = ykb
	   zh(nbod) = zkb
	   vxh(nbod) = vxkb
	   vyh(nbod) = vykb
	   vzh(nbod) = vzkb

	   vss = sqrt(vxkb**2 + vykb**2 + vzkb**2)
	   !bmin = RAT_T*vss*dt/365.0d0/TWOPI
	   !rdotv = vxkb*xkb+vykb*ykb+vzkb*zkb
	   !rmag = sqrt(xkb**2 + ykb**2 + zkb**2)
	   !costhet = rdotv/vss/rmag
	   !bmin = rmag * sqrt(1.d0-costhet*costhet)

	   r2hill(nbod) = 0.d0
	   rplsq(nbod) = (xkb*xkb+ykb*ykb+zkb*zkb) *1.01d0
       else
	   !just do impulse approximation if it's a small KBO
           call kb_enc_impapp(mkb,xkb,ykb,zkb,vxkb,vykb,
     &       vzkb,ntp,xht,yht,zht,vxht,vyht,vzht)
       endif
	   

       return
       end
