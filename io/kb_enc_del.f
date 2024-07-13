c------------------------------------------------------

       subroutine kb_enc_del(id,nbod,mass,xh,yh,zh,vxh,vyh,
     &           vzh,rplsq,icenc)

       include '../swift.inc'
       include '../stars2.inc'

       integer nbod,icenc,id,i
       real*8 mass(NPLMAX),rplsq(NPLMAX)
       real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
       real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
       real*8 r2hill(NPLMAX)

       common/c_rhill/r2hill

c-----
c...  Executable code 


       if(id.eq.nbod) then

         nbod = nbod - 1
         icenc = icenc - 1

       else

         do i=id,nbod-1
            xh(i) = xh(i+1)
            yh(i) = yh(i+1)
            zh(i) = zh(i+1)
            vxh(i) = vxh(i+1)
            vyh(i) = vyh(i+1)
            vzh(i) = vzh(i+1)
            mass(i) = mass(i+1)
            r2hill(i) = r2hill(i+1)
            rplsq(i) = rplsq(i+1)
         enddo
         nbod = nbod - 1
         icenc = icenc - 1

       endif

       return
       end
