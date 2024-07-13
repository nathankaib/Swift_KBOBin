import pylab
import struct
import glob
import math
import numpy as np
import os
from scipy import stats

def velgen(dvdist,npass):
    """this accepts a sample of encounter velocities and samples the CDF to return nenc velocities"""

    histdata = np.histogram(dvdist,bins=1000)
    histcum = np.cumsum(histdata[0])
    histcum = np.concatenate((np.asarray(0).flatten(),histcum))
    histcum = histcum/float(histcum[-1])
    histvals = np.concatenate((np.asarray(0).flatten(),histdata[1][1:]))

    #remove divisions where nothing changes
    deltacum = histcum[1:] - histcum[:-1]
    nochange = np.where(deltacum==0)[0]
    while (len(nochange)>0):
        histcum = np.delete(histcum,nochange+1)
        histvals = np.delete(histvals,nochange+1)
        deltacum = histcum[1:] - histcum[:-1]
        nochange = np.where(deltacum==0)[0]

    #now generate velocities from CDF sampling
    velvals = np.zeros(npass)
    toosmall = np.where(velvals < min(dvdist))[0]
    while (len(toosmall)>0):
        randnums = np.random.random_sample(len(toosmall))
        for i,histval in enumerate(histvals[:-1]):
            deltax = histvals[i+1] - histvals[i]
            deltay = histcum[i+1] - histcum[i]
            slope = deltay / deltax

            sel = np.where((randnums<=histcum[i+1])&(randnums>histcum[i]))[0]
            if len(sel) > 0:
                deltayval = randnums[sel] - histcum[i]
                velvals[toosmall[sel]] = histvals[i] + deltayval / slope

        toosmall = np.where(velvals < min(dvdist))[0]

    return(velvals)


def Nbod_calc(brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax):
    """this code calculates the number of bodies between two radii, Rvalidmin and Rvalidmax, assuming a continuous broken power-law SFD upto some maximum KBO size, Rbiggest, and an observed number of bodies above R=50 km, Npopabove50"""
    #inputs are the dirname:     directory name that contains encounter stats
    #               brightslope: the power-law index of the bright side of the luminosity distribution 
    #               faintslope:  the power-law index of the faint side of the luminosity distribution 
    #               rbreak:      the radius at which the bright-end SFD breaks to the faint-end SFD 
    #               Rbiggest:    the largest object within this particular KBO reservoir
    #               Npopabove50: the observationally inferred number of bodies in this particular reservoir with r>50 km
    #               Rvalidmin:   the minimum object radius we wish to generate
    #               Rvalidmax:   the maximum object radius we wish to generate

    #given magnitude slopes, calcluate radii slopes
    brightq = 5. * brightslope + 1
    faintq = 5. * faintslope + 1

    if rbreak >= 50: #have to consider both ends of size distribution then
        #first calculate the constant in front of the bright-end power-law
        brightk = Npopabove50 * (1. / (1. - faintq) * rbreak**(faintq - brightq) * (rbreak**(1. - faintq) - 50.**(1.-faintq)) + 1. / (1. - brightq) * (Rbiggest**(1. - brightq) - rbreak**(1. - brightq)))**(-1.)

        #now use requirement of continuity to get faint-end constant
        faintk = brightk * rbreak**(faintq - brightq)
    else: 
        #first calculate the constant in front of the bright-end power-law
        brightk = Npopabove50 * (1. / (1. - brightq) * (Rbiggest**(1. - brightq) - 50**(1. - brightq)))**(-1.)

        #now use requirement of continuity to get faint-end constant
        faintk = brightk * rbreak**(faintq - brightq)

    #now calculate the number of bodies within the user's specified size range
    if ((Rvalidmin < rbreak)&(Rvalidmax >= rbreak)&(Rvalidmax <= Rbiggest)):
        Nbod = faintk / (1. - faintq) * (rbreak**(1. - faintq) - Rvalidmin**(1. - faintq)) + brightk / (1. - brightq) * (Rvalidmax**(1. - brightq) - rbreak**(1. - brightq))

    if ((Rvalidmin < rbreak)&(Rvalidmax >= rbreak)&(Rvalidmax > Rbiggest)):
        Nbod = faintk / (1. - faintq) * (rbreak**(1. - faintq) - Rvalidmin**(1. - faintq)) + brightk / (1. - brightq) * (Rbiggest**(1. - brightq) - rbreak**(1. - brightq))

    if ((Rvalidmin < rbreak)&(Rvalidmax < rbreak)):
        Nbod = faintk / (1. - faintq) * (Rvalidmax**(1. - faintq) - Rvalidmin**(1. - faintq))

    if ((Rvalidmin >= rbreak)&(Rvalidmax >= rbreak)&(Rvalidmax <= Rbiggest)):
        Nbod = brightk / (1. - brightq) * (Rvalidmax**(1. - brightq) - Rvalidmin**(1. - brightq))

    if ((Rvalidmin >= rbreak)&(Rvalidmax >= rbreak)&(Rvalidmax > Rbiggest)):
        Nbod = brightk / (1. - brightq) * (Rbiggest**(1. - brightq) - Rvalidmin**(1. - brightq))

    return(Nbod,brightk,faintk)

def R_generator(randvals, brightslope, faintslope, rbreak, brightk, faintk, Rvalidmin, Rvalidmax, Nbod):
    """given slopes, break, and SFD constants, this code generates random KBO radii from Rvalidmin to Rvalidmax assuming there are Nbod bodies if you integrate SFD from Rvalidmin to Rvalidmax"""
    #inputs are the randvals:    random seeds for function that inverts the CDF of SFD
    #               brightslope: the power-law index of the bright side of the luminosity distribution 
    #               faintslope:  the power-law index of the faint side of the luminosity distribution 
    #               rbreak:      the radius at which the bright-end SFD breaks to the faint-end SFD 
    #               brightk:     the bright-end SFD coefficient
    #               faintk:      the faint-end SFD coefficient
    #               Rvalidmin:   the minimum object radius we wish to generate
    #               Rvalidmax:   the maximum object radius we wish to generate
    #               Nbod:        the number of KBOs between Rvalidmin and Rvalidmax

    #Note: these 2 parameters (encounter distance and number of bodies) can be changed when generating encounters.txt files for a population(s)
    #but so far, I've kept them at 10,000 and .01 AU

    #given magnitude slopes, calcluate radii slopes
    brightq = 5. * brightslope + 1
    faintq = 5. * faintslope + 1

    radvals = np.zeros(len(randvals))
    if ((Rvalidmin<=rbreak)&(Rvalidmax>=rbreak)):
        integraltobreak = faintk / (Nbod * (1. - faintq)) * (rbreak**(1. - faintq) - Rvalidmin**(1. - faintq))
        sel = np.where(randvals < integraltobreak)[0]
        radvals[sel] = (rbreak**(1. - faintq) - randvals[sel] * Nbod * (1. - faintq) / faintk)**(1. / (1. - faintq))
        sel = np.where(randvals >= integraltobreak)[0]
        radvals[sel] = (Nbod * (1. - brightq) / brightk * (randvals[sel] - faintk / Nbod / (1. - faintq) * (rbreak**(1. - faintq) - Rvalidmin**(1. - faintq))) + rbreak**(1. - brightq))**(1. / (1. - brightq))

    if ((Rvalidmin<=rbreak)&(Rvalidmax<=rbreak)):
        radvals = (randvals * Nbod * (1. - faintq) / faintk + Rvalidmin**(1. - faintq))**(1. / (1. - faintq))

    if (Rvalidmin>=rbreak):
        radvals = (randvals * Nbod * (1. - brightq) / brightk + Rvalidmin**(1. - brightq))**(1. / (1. - brightq))

    return(radvals)


def passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength):
    """this routine generates the info for KBO passages from a KBO population"""
    #inputs are the dirname:     directory name that contains encounter stats
    #               brightslope: the power-law index of the bright side of the luminosity distribution 
    #               faintslope:  the power-law index of the faint side of the luminosity distribution 
    #               rbreak:      the radius at which the bright-end SFD breaks to the faint-end SFD 
    #               Rbiggest:    the largest object within this particular KBO reservoir
    #               Npopabove50: the observationally inferred number of bodies in this particular reservoir with r>50 km
    #               Rvalidmin:   the minimum object radius we wish to generate
    #               Rvalidmax:   the maximum object radius we wish to generate
    #               RSphere:     the radius at which KBO passages will begin
    #               tlength:     the amount of time for which we will generate KBO passages

    #Note: these 2 parameters (encounter distance and number of bodies) can be changed when generating encounters.txt files for a population(s)
    #but so far, I've kept them at 10,000 and .01 AU
    Nparts = 1e4
    collstatrad = 2.5e-3 * 4. * 1.5e8 #in km

    #reading encounter sim info
    finishedsimlist = glob.glob(dirname+'/*/finished')
    if dirname == '../CCB_CCB':
        Npairs = (Nparts**2./2. + Nparts/2.) * len(finishedsimlist)
    else:
        Npairs = Nparts**2. * len(finishedsimlist)

    #read encounter values
    encfilename = dirname + '/encounters_all.txt'
    a1,ecc1,inc1,a2,ecc2,inc2,dv,smin = np.genfromtxt(encfilename,unpack=True,skip_header=2)

    #get encounter times
    nencs = len(dv)
    coprob = nencs / 323. / Npairs / collstatrad**2. #in cols/yr/km^2

    #now get the total number of bodies between Rvalidmin and Rvalidmax
    #as well as the SFD coefficients (the k in dN/dr = k * r^-q)
    Nbod,brightk,faintk = Nbod_calc(brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax)

    #use RSphere, Nbod, collision prob and total time length to get # of KBO passages
    npass = coprob * Nbod * RSphere**2. * tlength
    leftover = np.mod(npass,int(npass))
    testforextra = np.random.random_sample(1)
    if testforextra<=leftover:
        npass= int(npass)+1
    else:
        npass= int(npass)

    #passage times
    times = np.random.random_sample(npass) * tlength

    #now generate seeds for radius generating routine
    randvals = np.random.random_sample(npass)

    #use seeds to generate radii
    radvals = R_generator(randvals, brightslope, faintslope, rbreak, brightk, faintk, Rvalidmin, Rvalidmax, Nbod)

    #get encounter velocities
    vmagall = velgen(dv,npass)

    #set initial distances of passing objects
    rstartall = RSphere + vmagall * 0.

    #returns passage times, KBO sizes, start distances, and velocity magnitudes
    return(times,radvals,rstartall,vmagall)

def adjust_velocity_directions(vmagall):
    """this function generates random velocity vectors for KBO encounters penetrating an imaginary sphere around a binary"""
    #here we make sure that the ratio of on-axis vs grazing encounters is right (see Henon 1972)

    #we're implicitly assuming that all passages start at the positive z-axis
    #we'll randomize position and velocity vectors across the imaginary sphere later on

    #computing direction stars are moving in when they enter sphere.
    #the z-axis for the coordinate systems of these directions are perpendicular    
    #to the inner surface of the sphere at the star's point of entry (Henon 1972)
    #(confining all velocities to a 2-D plane initially)
    npass=len(vmagall)
    vincu=np.random.random_sample(npass)
    vinc=np.arcsin(np.sqrt(vincu))

    #now i'm adding a random positive or negative sign to the velocity angle to make
    #the distribution run from -pi/2 to pi/2
    sign=np.random.random_sample(npass)
    posi=np.unique(np.asarray(np.where(sign > .5)))
    neg=np.unique(np.asarray(np.where(sign < .5)))

    sign[posi] = 1.0
    sign[neg] = -1.0

    vinc=vinc*sign

    #all motion is confined to the x-z plane initially
    vx=vmagall*np.sin(vinc)
    vz=-vmagall*np.cos(vinc)
    vy=vz*0.

    return(vx, vy, vz)

def randomize_posvel(rstartall, vx, vy, vz):
    """this function takes position and velocity vectors of KBO passages and randomly distributes them across sphere surrounding binary"""
    #here we assume that all passages are starting at postive z-axis

    npass=len(rstartall)
    #see!
    #initially having every star enter sphere at same point
    x, y, z = np.zeros(npass),np.zeros(npass),np.zeros(npass)
    x[0:npass]=0.
    y[0:npass]=0.
    z[0:npass]=rstartall

    #rotating the point-of-entry position and velocity vectors around the
    #z-axis by a random angle to create a 3-dimensional distribution of impact
    #angles relative to the point of entry (a rotated velocity vector
    #would trace out a cone overhead)

    pos = np.zeros((3,1,npass))
    vel = np.zeros((3,1,npass))

    pos[0,0,:] = x+0.0
    pos[1,0,:] = y+0.0
    pos[2,0,:] = z+0.0
    vel[0,0,:] = vx+0.0
    vel[1,0,:] = vy+0.0
    vel[2,0,:] = vz+0.0

    phi=np.random.random_sample(npass)*2.*math.pi
    rotz = np.zeros((3,3,npass))
    rotz[0,0,:] = np.cos(phi)
    rotz[0,1,:] = -np.sin(phi)
    rotz[0,2,:] = 0.
    rotz[1,0,:] = np.sin(phi)
    rotz[1,1,:] = np.cos(phi)
    rotz[1,2,:] = 0.
    rotz[2,0,:] = 0.
    rotz[2,1,:] = 0.
    rotz[2,2,:] = 1.

    pos = (rotz[:,:,None]*pos).sum(axis=1)
    vel = (rotz[:,:,None]*vel).sum(axis=1)

    #rotating point-of-entry position and velocity vectors around the x-axis by 
    #a random angle (here we make a half-ring of points-of-entry around the sun 
    #perpindicular to the x-axis)
    cosphi=np.random.random_sample(npass)*2.-1.
    phi=np.arccos(cosphi)

    rotx = np.zeros((3,3,npass))
    rotx[0,0,:] = 1.
    rotx[0,1,:] = 0.
    rotx[0,2,:] = 0.
    rotx[1,0,:] = 0.
    rotx[1,1,:] = np.cos(phi)
    rotx[1,2,:] = -np.sin(phi)
    rotx[2,0,:] = 0.
    rotx[2,1,:] = np.sin(phi)
    rotx[2,2,:] = np.cos(phi)

    pos = (rotx[:,:,None]*pos).sum(axis=1)
    vel = (rotx[:,:,None]*vel).sum(axis=1)

    #rotating point-of-entry position and velocity vectors around the z-axis by a
    #a random angle again (for each point-of-entry, this will twist our ring by a 
    #randomly selected angle between 0 and 2pi, which will finally generate a 3-d 
    #sphere of stellar passages in an isotropic field of stars)
    phi=np.random.random_sample(npass)*2.*math.pi
    rotz = np.zeros((3,3,npass))

    rotz[0,0,:] = np.cos(phi)
    rotz[0,1,:] = -np.sin(phi)
    rotz[0,2,:] = 0.
    rotz[1,0,:] = np.sin(phi)
    rotz[1,1,:] = np.cos(phi)
    rotz[1,2,:] = 0.
    rotz[2,0,:] = 0.
    rotz[2,1,:] = 0.
    rotz[2,2,:] = 1.

    pos = (rotz[:,:,None]*pos).sum(axis=1)
    vel = (rotz[:,:,None]*vel).sum(axis=1)

    #done with rotations
    x = pos[0,0,:]+0.0
    y = pos[1,0,:]+0.0
    z = pos[2,0,:]+0.0
    vx = vel[0,0,:]+0.0
    vy = vel[1,0,:]+0.0
    vz = vel[2,0,:]+0.0

    #transforming positions from km to AU 
    #velocities are already in AU/day
    x=x/1.50e8
    y=y/1.50e8
    z=z/1.50e8

    return(x, y, z, vx, vy, vz)


def passgen_all_new(passname,tlength):
    density = .5 #g/cm^3
    Msys = 1.0e21 #in g
    Msun = 1.9891e33 #in g
    asys = 44. #in AU
    QW322Hill = (2e21 / Msun / 3.)**.3333333 * asys * 1.5e8 #in km
    PlutoHill = (1.3e25 / Msun / 3.)**.33333333 * asys * 1.5e8 #in km
    #rbreak = 85. #break for Lawler knee
    rbreak = 65. #break for Lawler divot
    #rbreak = 50. #break for Lawler divot

    timeall, radall, rstartall, vmagall, impapproxall = np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0, dtype='int')

    #first do cold classical belt
    dirname = '../CCB_CCB'
    brightslope = 1.2
    faintslope = 0.4
    Rbiggest = 500.
    Npopabove50 = 95000.
    #do CCBs between 50 km and 150 km
    RSphere = 3.*QW322Hill #only look at within 3 hill radii, since there are many of them
    Rvalidmin = 50.
    Rvalidmax = 150.
    CCBtimes, CCBradvals, CCBrstarts, CCBvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    CCBimpapproxall = np.ones(len(CCBvmags))
    timeall = np.concatenate((timeall, CCBtimes))
    radall = np.concatenate((radall, CCBradvals))
    rstartall = np.concatenate((rstartall, CCBrstarts))
    vmagall = np.concatenate((vmagall, CCBvmags))
    impapproxall = np.concatenate((impapproxall,np.ones(len(CCBrstarts),dtype='int').flatten()))
    #do CCBs between 5 km and 50 km
    RSphere = QW322Hill #only look at within 3 hill radius, since there are many, many of them
    Rvalidmin = 5.
    Rvalidmax = 50.
    CCBtimes, CCBradvals, CCBrstarts, CCBvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    CCBimpapproxall = np.ones(len(CCBvmags))
    timeall = np.concatenate((timeall, CCBtimes))
    radall = np.concatenate((radall, CCBradvals))
    rstartall = np.concatenate((rstartall, CCBrstarts))
    vmagall = np.concatenate((vmagall, CCBvmags))
    impapproxall = np.concatenate((impapproxall,np.ones(len(CCBrstarts),dtype='int').flatten()))
    #do CCBs between 150 km and 500 km
    RSphere = QW322Hill #only look at within 3 hill radius, since there are many, many of them
    RSphere = 3.*PlutoHill #look at within 3 Pluto Hill radii, since they are big, and they are rare
    Rvalidmin = 150.
    Rvalidmax = 500.
    CCBtimes, CCBradvals, CCBrstarts, CCBvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    CCBimpapproxall = np.ones(len(CCBvmags))
    timeall = np.concatenate((timeall, CCBtimes))
    radall = np.concatenate((radall, CCBradvals))
    rstartall = np.concatenate((rstartall, CCBrstarts))
    vmagall = np.concatenate((vmagall, CCBvmags))
    impapproxall = np.concatenate((impapproxall,np.zeros(len(CCBrstarts),dtype='int').flatten())) #have to actually integrate these since encounter timescale is long

    #next do hot classical belt
    dirname = '../CCB_HCB'
    brightslope = 0.9
    faintslope = 0.4
    Rbiggest = 800.
    Npopabove50 = 35000.
    #do HCBs between 50 km and 150 km
    RSphere = 3.*QW322Hill #only look at within 3 hill radii, since there are many of them
    Rvalidmin = 50.
    Rvalidmax = 150.
    HCBtimes, HCBradvals, HCBrstarts, HCBvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    HCBimpapproxall = np.ones(len(HCBvmags))
    timeall = np.concatenate((timeall, HCBtimes))
    radall = np.concatenate((radall, HCBradvals))
    rstartall = np.concatenate((rstartall, HCBrstarts))
    vmagall = np.concatenate((vmagall, HCBvmags))
    impapproxall = np.concatenate((impapproxall,np.ones(len(HCBrstarts),dtype='int').flatten()))
    #do HCBs between 5 km and 50 km
    #RSphere = QW322Hill #only look at within 3 hill radius, since there are many, many of them
    #Rvalidmin = 5.
    #Rvalidmax = 50.
    #HCBtimes, HCBradvals, HCBrstarts, HCBvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    #HCBimpapproxall = np.ones(len(HCBvmags))
    #timeall = np.concatenate((timeall, HCBtimes))
    #radall = np.concatenate((radall, HCBradvals))
    #rstartall = np.concatenate((rstartall, HCBrstarts))
    #vmagall = np.concatenate((vmagall, HCBvmags))
    #impapproxall = np.concatenate((impapproxall,np.ones(len(HCBrstarts),dtype='int').flatten()))
    #do HCBs between 150 km and 800 km
    RSphere = QW322Hill #only look at within 3 hill radius, since there are many, many of them
    RSphere = 3.*PlutoHill #look at within 3 Pluto Hill radii, since they are big, and they are rare
    Rvalidmin = 150.
    Rvalidmax = 800.
    HCBtimes, HCBradvals, HCBrstarts, HCBvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    HCBimpapproxall = np.ones(len(HCBvmags))
    timeall = np.concatenate((timeall, HCBtimes))
    radall = np.concatenate((radall, HCBradvals))
    rstartall = np.concatenate((rstartall, HCBrstarts))
    vmagall = np.concatenate((vmagall, HCBvmags))
    impapproxall = np.concatenate((impapproxall,np.zeros(len(HCBrstarts),dtype='int').flatten())) #have to actually integrate these since encounter timescale is long

    #next do Plutino population
    dirname = '../CCB_3to2'
    brightslope = 0.9
    faintslope = 0.4
    Rbiggest = 1200.
    Npopabove50 = 13000.
    #do Plutinos between 50 km and 150 km
    RSphere = 3.*QW322Hill #only look at within 3 hill radii, since there are many of them
    Rvalidmin = 50.
    Rvalidmax = 150.
    Pluttimes, Plutradvals, Plutrstarts, Plutvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    Plutimpapproxall = np.ones(len(Plutvmags))
    timeall = np.concatenate((timeall, Pluttimes))
    radall = np.concatenate((radall, Plutradvals))
    rstartall = np.concatenate((rstartall, Plutrstarts))
    vmagall = np.concatenate((vmagall, Plutvmags))
    impapproxall = np.concatenate((impapproxall,np.ones(len(Plutrstarts),dtype='int').flatten()))
    #do Plutinos between 5 km and 50 km
    #RSphere = QW322Hill #only look at within 3 hill radius, since there are many, many of them
    #Rvalidmin = 5.
    #Rvalidmax = 50.
    #Pluttimes, Plutradvals, Plutrstarts, Plutvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    #Plutimpapproxall = np.ones(len(Plutvmags))
    #timeall = np.concatenate((timeall, Pluttimes))
    #radall = np.concatenate((radall, Plutradvals))
    #rstartall = np.concatenate((rstartall, Plutrstarts))
    #vmagall = np.concatenate((vmagall, Plutvmags))
    #impapproxall = np.concatenate((impapproxall,np.ones(len(Plutrstarts),dtype='int').flatten()))
    #do Plutinos between 150 km and 1200 km
    RSphere = QW322Hill #only look at within 3 hill radius, since there are many, many of them
    RSphere = 3.*PlutoHill #look at within 3 Pluto Hill radii, since they are big, and they are rare
    Rvalidmin = 150.
    Rvalidmax = 1200.
    Pluttimes, Plutradvals, Plutrstarts, Plutvmags = passagepop_gen(dirname, brightslope, faintslope, rbreak, Rbiggest, Npopabove50, Rvalidmin, Rvalidmax, RSphere, tlength)
    Plutimpapproxall = np.ones(len(Plutvmags))
    timeall = np.concatenate((timeall, Pluttimes))
    radall = np.concatenate((radall, Plutradvals))
    rstartall = np.concatenate((rstartall, Plutrstarts))
    vmagall = np.concatenate((vmagall, Plutvmags))
    impapproxall = np.concatenate((impapproxall,np.zeros(len(Plutrstarts),dtype='int').flatten())) #have to actually integrate these since encounter timescale is long

    #now order everything by passage time
    order = np.argsort(timeall)
    timeall = timeall[order]
    radall = radall[order]
    rstartall = rstartall[order]
    vmagall = vmagall[order]
    impapproxall = impapproxall[order]
    
    #now orient velocity vectors to attain the right impact parameter distribution
    vx, vy, vz = adjust_velocity_directions(vmagall)

    #now orient position and velocity vectors across an imaginary sphere surrounding binary 
    x, y, z, vx, vy, vz = randomize_posvel(rstartall, vx, vy, vz)

    #now use density and radii to calculate masses
    randmassall = 4./3. * math.pi * density * (radall*1e5)**3. #in g
    print(max(np.isnan(randmassall)))

    vmag = np.sqrt(vx**2.+vy**2. + vz**2.)
    tbar = ((0.-x)*vx + (0.-y)*vy + (0.-z)*vz) / (vmag*vmag)
    xclosebar = x + vx*tbar
    yclosebar = y + vy*tbar
    zclosebar = z + vz*tbar
    rbar = np.sqrt(xclosebar**2. + yclosebar**2. + zclosebar**2.)
    tenc = rbar / vmag

    f = open(passname,'wb')
    print('start write')
    oneline = struct.pack('2hi2h',*tuple((int(4),0,int(len(timeall)),int(4),0)))
    f.write(oneline)
    for i in range(len(timeall)):
        oneline = struct.pack('2h9f2h',*tuple((int(4+8*4),0,tenc[i],timeall[i]*365.25,randmassall[i]/1.9891e33*2.95913978e-4,x[i],y[i],z[i],vx[i],vy[i],vz[i],int(4+8*4),0)))
        f.write(oneline)
    f.close()

    return()


for i in range(1,2):
    print(i)
    filename = 'passages' + str(i) + '.in'
    passgen_all_new(filename,1e8)

