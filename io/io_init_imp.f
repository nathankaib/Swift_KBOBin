c---------------------------------------------------------------------
        subroutine io_init_imp(impfile,nimp,iimp,timp,impobj,
     &     mimp,ximp,yimp,zimp,vximp,vyimp,vzimp,tstop,t0)


       include '../swift.inc'
       include 'io.inc'
       include '../stars2.inc'

c...  Input
       character(*) impfile
       real*8 tstop,t0

c...  Output
       real*8 timp(NTPMAX),ximp(NTPMAX),yimp(NTPMAX),zimp(NTPMAX)
       real*8 vximp(NTPMAX),vyimp(NTPMAX),vzimp(NTPMAX),mimp(NTPMAX)
       integer iimp,nimp,iimpt,impobj(NTPMAX)

c...   Internal
       integer i,ierr,test
       real*8 itype

c-----
c...  Executable code      

       write(*,*) 'KBO impact data file is ',impfile
       call io_open(7,impfile,'old','formatted',ierr)
       if(ierr.ne.0) then
         write(*,*) 'Could not open KBO encounter data file'
         call util_exit(1)
       endif

       read(7,*) nimp
       write(*,*) 'There are ',nimp,' KBO impacts'
       do i=1,nimp

c     old format
         read(7,*) timp(i),impobj(i),mimp(i),ximp(i),yimp(i),
     %       zimp(i),vximp(i),vyimp(i),vzimp(i)
         if(timp(i).gt.tstop) then
            goto 1
         endif
       enddo
       
c     if we use all the encounters
       if (i.gt.nimp) then
          timp(i) = 1d80 !set next absurdly high so we never use it
       endif

 1     nimp = i
       write(*,*) '    Will use ',nimp,' of them.'

       if(timp(nimp).lt.tstop) then
         write(*,*) 'KBO impact file is too short'
         write(*,*) 'tstop = ',tstop,' timp(nimp) = ',timp(nimp)
         call util_exit(1)
       endif

       iimp = -1
       do iimpt=1,nimp
         if((timp(iimpt).gt.t0).and.(iimp.lt.0)) then
            iimp = iimpt
         endif
       enddo

       write(*,*) 'The integration starts at',t0
       write(*,*) 'The first impact will be at ',timp(iimp)
       if(iimp.eq.1) then
         write(*,*) iimp
       else
         write(*,*) iimp,timp(iimp),timp(iimp-1)
       endif
       return
       end
