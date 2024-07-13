c-----------------------------------------------------

       subroutine kb_enc_chk(nbod,nbod0,mass,xh,yh,zh,vxh,vyh,
     &     vzh,rplsq,icenc,i1st,rmax)

       include '../swift.inc'
       include '../stars2.inc'

       integer nbod,icenc,i,i1st,nbod0
       real*8 mass(NPLMAX),r2max,r2,rplsq(NPLMAX),rmax
       real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
       real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c-----
c...  Executable code 

       do i = nbod0+1, nbod

         r2 = (xh(i)**2 + yh(i)**2 + zh(i)**2)

         if(r2.ge.rplsq(i)) then
            call kb_enc_del(i,nbod,mass,xh,yh,zh,vxh,vyh,
     &           vzh,rplsq,icenc)
	    write(*,*) 'removing body'
	    write(*,*) mass(i)
            i1st = 0
         endif

       enddo

       return
       end
