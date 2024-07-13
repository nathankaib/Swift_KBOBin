       subroutine kb_enc_impapp(mkb,xkb,ykb,zkb,vxkb,vykb,
     &       vzkb,ntp,xht,yht,zht,vxht,vyht,vzht)


       include '../swift.inc'
       include '../stars2.inc'

       integer ntp,i
       real*8 xkb,ykb,zkb,vxkb,vykb,vzkb,mkb
       real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
       real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

       real*8 vmag,t,xclose,yclose,zclose
       real*8 dx,dy,dz,dist,dunitx,dunity,dunitz
       real*8 dsun,dsununitx,dsununity,dsununitz
       real*8 impfactor,impsunx,impsuny,impsunz
       real*8 impx,impy,impz

c-----
c...  Executable code 
       vmag = sqrt(vxkb*vxkb + vykb*vykb + vzkb*vzkb)
       do i = 1,ntp
          !figuring out time of closest approach to binary companion
	  t = ((xht(i)-xkb)*vxkb + (yht(i)-ykb)*vykb + 
     &         (zht(i)-zkb)*vzkb) / (vmag*vmag)
	  xclose = xkb + vxkb*t
	  yclose = ykb + vykb*t
	  zclose = zkb + vzkb*t

	  !figuring out actual closest approach distance to binary
	  dx = xclose - xht(i)
	  dy = yclose - yht(i)
	  dz = zclose - zht(i)

	  !figuring out unit vector from companion to passing star
	  dist = sqrt(dx*dx + dy*dy + dz*dz)
	  dunitx = dx/dist
	  dunity = dy/dist
	  dunitz = dz/dist

	  !figuring out time of closest approach to sun (reusing t var.)
	  t = (-1.0)*(xkb*vxkb + ykb*vykb + zkb*vzkb) / (vmag*vmag)

	  !figuring out actual closest approach to sun (reusing vars.)
	  dx = xkb + vxkb*t
	  dy = ykb + vykb*t
	  dz = zkb + vzkb*t

	  !figuring out unit vector from sun to passing star
	  dsun = sqrt(dx*dx + dy*dy + dz*dz)
	  dsununitx = dx/dsun
	  dsununity = dy/dsun
	  dsununitz = dz/dsun

	  !calculating impulse on sun
	  impfactor = 2.d0 * mkb / vmag
	  impsunx = impfactor*dsununitx/dsun
	  impsuny = impfactor*dsununity/dsun
	  impsunz = impfactor*dsununitz/dsun
	   
	  !calculating relative impulse on binary companion
	  impx = impfactor*dunitx/dist - impsunx
	  impy = impfactor*dunity/dist - impsuny
	  impz = impfactor*dunitz/dist - impsunz
	   
	  !adding impulse to binary's velocity
	  vxht(i) = vxht(i) + impx
	  vyht(i) = vyht(i) + impy
	  vzht(i) = vzht(i) + impz
       enddo

       return
       end
