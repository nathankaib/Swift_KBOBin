c**********************************************************************
c		      SWIFT_RMVS4.F
c**********************************************************************
c
c                 INCLUDES CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c
c Authors: Nate Kaib 
c Date:    1/28/19
c Last revision: 1/28/19

     
	include 'swift.inc'
	include '../stars2.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT),i1st
	integer nbod,ntp,nleft,npl
	integer iflgchk,iub,iuj,iud,iue
        real*8 rstat(NTPMAX,NSTATR)

	real*8 tkb(NSTENC),xkb(NSTENC),ykb(NSTENC),zkb(NSTENC)
	real*8 vxkb(NSTENC),vykb(NSTENC),vzkb(NSTENC),mkb(NSTENC)
        real*8 tenckb(NSTENC)
        integer impappflg

	integer ikb,nbod0,icenc,nkb
	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat
	character*280 kbfile


c-----
c...    Executable code 

c...    print version number
        call util_version

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,npl,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

c Get data for the KBO passages
	write(*,*) 'Enter name of KBO encounter file : '
        read(*,999) kbfile
        call io_init_kb_bin(kbfile,nkb,ikb,tenckb,tkb,mkb,xkb,ykb,
     &     zkb,vxkb,vykb,vzkb,tstop,t0)

c Initialize solar orbit settings
	call accsun_init

c Initialize initial time and times for first output and first dump
	t = t0
	tout = t0 + dtout
	tdump = t0 + dtdump

        iub = 20
        iuj = 30
        iud = 40
        iue = 60

c...    Do the initial io write
        if(btest(iflgchk,0))  then ! bit 0 is set
           call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,1))  then ! bit 1 is set
           call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,2))  then    ! bit 2 is set
           eoff = 0.0d0
           call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &          vyh,vzh,iue,fopenstat,eoff)
        endif
        if(btest(iflgchk,3))  then    ! bit 3 is set
           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
        endif

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           call io_discard_write(0,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &          vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &          'discard.out',fopenstat,nleft)
        endif

        nleft = ntp
	nbod0=npl
	i1st = 0
c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '
	  do while ( (t .le. tstop) .and. 
     &       ((ntp.eq.0).or.(nleft.gt.0)) )

	     if (nbod.gt.1) then
	        write(*,*) xh(2),yh(2),zh(2)
	     endif


 	     call rmvs4_step(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &         xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &         vzht,istat,rstat,dt)
	     if (nbod.gt.1) then
	        write(*,*) xh(2),yh(2),zh(2)
	     endif

	     t = t + dt

             if(btest(iflgchk,4))  then    ! bit 4 is set
                call discard(t,dt,nbod,npl,ntp,mass,xh,yh,zh,vxh,vyh,
     &               vzh,xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
     &               qmin,lclose,rplsq,istat,rstat)
                call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &               vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &               'discard.out',fopenstat,nleft)
             else
                nleft = ntp
             endif

	     if(icenc.ne.0) then
		!note: rplsq is now used to decide when to remove KBO passers
		call kb_enc_chk(nbod,nbod0,mass,xh,yh,zh,vxh,vyh,vzh,
     &            rplsq,icenc,i1st,rmax)
	     endif
c start a stellar encounter if it is time
	     if(t.gt.tkb(ikb)) then
		if (tenckb(ikb).gt.dt) then
		    impappflg = 0 !long encounter so do direct integration
		else
		    impappflg = 1 !short encounter so do impulse approx.
		endif

		call kb_enc_add(nbod,mass,xh,yh,zh,vxh,vyh,vzh,rplsq,
     &           mkb(ikb),xkb(ikb),ykb(ikb),zkb(ikb),vxkb(ikb),
     &           vykb(ikb),vzkb(ikb),dt,impappflg,icenc,ntp,xht,
     &           yht,zht,vxht,yht,vxht,vyht,vzht)
		ikb = ikb + 1
		i1st = 0
		write(*,*) 'start encounter'
	     endif



c if it is time, output orb. elements, 
	  if(t .ge. tout) then 

             if(btest(iflgchk,0))  then    ! bit 0 is set
                call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &               iub,fopenstat)
              endif
             if(btest(iflgchk,1))  then    ! bit 1 is set
                call  io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &               iub,fopenstat)
              endif

	    tout = tout + dtout
	  endif

c If it is time, do a dump
          if(t.ge.tdump) then

             tfrac = (t-t0)/(tstop-t0)
             write(*,998) t,tfrac,nleft
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of active tp =',i4)
	     call io_dump_pl('dump_pl.dat',nbod,npl,mass,xh,yh,zh,
     &                 vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	     call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &                      vxht,vyht,vzht,istat,rstat)
	     call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
	     tdump = tdump + dtdump

             if(btest(iflgchk,2))  then    ! bit 2 is set
                call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &               xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
             endif
             if(btest(iflgchk,3))  then    ! bit 3 is set
                call anal_jacobi_write(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,
     &               iuj,fopenstat)
             endif

	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

        call io_dump_pl('dump_pl.dat',nbod,npl,mass,xh,yh,zh,
     &            vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &              vxht,vyht,vzht,istat,rstat)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_rmvs4.f
c---------------------------------------------------------------------
