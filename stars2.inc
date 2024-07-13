c*************************************************************************
c                        STARS2.INC
c*************************************************************************
C Include file for the stellar encounter stuff.
c
c Author:  Hal Levison  
c Date:    10/16/98
c Last revision: 

c...   Maximum number of stellar encounters
       integer  NSTENC
       parameter  (NSTENC = 20000000) 

c...   Ratio of timestep to `orbital period' to resolve peri
       real*8 RAT_T
       parameter (RAT_T=600.0d0)  

