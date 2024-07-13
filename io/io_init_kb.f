c---------------------------------------------------------------------
        subroutine io_init_kb(kbfile,nkb,ikb,tenckb,tkb,mkb,xkb,ykb,
     &               zkb,vxkb,vykb,vzkb,tstop,t0)

       include '../swift.inc'
       include 'io.inc'
       include '../stars2.inc'

c...  Input
       character(*) kbfile
       real*8 tstop,t0

c...  Output
       real*8 tkb(NSTENC),xkb(NSTENC),ykb(NSTENC),zkb(NSTENC)
       real*8 vxkb(NSTENC),vykb(NSTENC),vzkb(NSTENC),mkb(NSTENC)
       real*8 tenckb(NSTENC)
       integer ikb,nkb,ikbt

c...   Internal
       integer i,ierr,test
       real*8 itype

c-----
c...  Executable code      

       write(*,*) 'KBO encounter data file is ',kbfile
       call io_open(7,kbfile,'old','formatted',ierr)
       if(ierr.ne.0) then
         write(*,*) 'Could not open KBO encounter data file'
         call util_exit(1)
       endif

       read(7,*) nkb
       write(*,*) 'There are ',nkb,' KBO encounters'
       do i=1,nkb

c     old format
         read(7,*) tenckb(i),tkb(i),mkb(i),xkb(i),ykb(i),
     %       zkb(i),vxkb(i),vykb(i),vzkb(i)
         if(tkb(i).gt.tstop) then
            goto 1
         endif
       enddo
       
c     if we use all the encounters
       if (i.gt.nkb) then
          tkb(i) = 1d80 !set next absurdly high so we never use it
       endif

 1     nkb = i
       write(*,*) '    Will use ',nkb,' of them.'

       if(tkb(nkb).lt.tstop) then
         write(*,*) 'KBO encounter file is too short'
         write(*,*) 'tstop = ',tstop,' tkb(nkb) = ',tkb(nkb)
         call util_exit(1)
       endif

       ikb = -1
       do ikbt=1,nkb
         if((tkb(ikbt).gt.t0).and.(ikb.lt.0)) then
            ikb = ikbt
         endif
       enddo

       write(*,*) 'The integration starts at',t0
       write(*,*) 'The first encounter will be at ',tkb(ikb)
       if(ikb.eq.1) then
         write(*,*) ikb
       else
         write(*,*) ikb,tkb(ikb),tkb(ikb-1)
       endif
       return
       end
