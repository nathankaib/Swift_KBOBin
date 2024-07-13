c---------------------------------------------------------------------
c       Subroutine to setup the galactic tides ACCSUN_INIT
c---------------------------------------------------------------------
c This subroutine reads in the constants for the Sun's acceleration on
c the system and puts them in a common block
c 
c             Input:  NONE
c
c             Output: NONE
c
c             Common:
c  msun                  ==>  mass of sun
c
c  ahel                  ==>  semimajor axis of orbit about Sun in AU
c
c NOTE: !!!!!!!!!!!!!
c                    
c
c Authors:  Nathan Kaib
c Date:    2/07/18
c Last revision: 2/07/18

       subroutine accsun_init

       include '../../swift.inc'

c...  Common
       real*8 msun, ahel

c...  Internal
       common/c_sunacc/msun,ahel

c----
c...  Executable code 

       write(*,*) ' Mass of central star: '
       read(*,*) msun
       write(*,*) msun
       write(*,*) ' Heliocentric semimajor axis (in AU): '
       read(*,*) ahel
       write(*,*) ahel

       return
       end
