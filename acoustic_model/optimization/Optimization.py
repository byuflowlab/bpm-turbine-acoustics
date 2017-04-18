from pyoptsparse import Optimization, SNOPT
import numpy as np
from numpy import fabs, pi, sin, cos
from scipy.optimize import fsolve
from ccblade import CCAirfoil, CCBlade
from scipy.interpolate import Akima1DInterpolator
from os import path
from sys import argv

import _florisunify
import _bpmacoustic

from joblib import Parallel, delayed

# PARAMETERS USED FOR TUNING THE FLORIS MODEL
def floris_parameters():
    pP = 1.88
    ke = 0.065
    keCorrArray = 0.0
    keCorrCT = 0.0
    Region2CT = 4.0*(1.0/3.0)*(1.0-(1.0/3.0))
    kd = 0.15
    # me = np.array([-0.5, 0.22, 1.0]) # cosine off
    me = np.array([-0.5, 0.3, 1.0]) # cosine on

    initialWakeDisplacement = -4.5
    initialWakeAngle = 1.5

    baselineCT = 4./3.*(1.-1./3.)

    keCorrTI = 0.0
    baselineTI = 0.045

    keCorrHR = 0.0 # neutral, with heating rate 0, is baseline
    keCorrHRTI = 0.0
    keSaturation = 0.0

    kdCorrYawDirection = 0.0

    MU = np.array([0.5, 1.0, 5.5])

    CTcorrected = True # CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)
    CPcorrected = True # CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)

    axialIndProvided = True # CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)

    # useWakeAngle = False # cosine off
    useWakeAngle = True # cosine on

    bd = -0.01

    useaUbU = True
    aU = 5.0 # degrees
    bU = 1.66

    adjustInitialWakeDiamToYaw = False

    FLORISoriginal = False # override all parameters and use FLORIS as original in first Wind Energy paper

    # cos_spread = 1E12 # spread of cosine smoothing factor (percent of sum of wake and rotor radii) cosine off
    cos_spread = 2.0 # spread of cosine smoothing factor (percent of sum of wake and rotor radii) cosine on

    return pP, ke, keCorrArray, keCorrCT, Region2CT, kd, me, initialWakeDisplacement, initialWakeAngle, baselineCT, keCorrTI, baselineTI, keCorrHR, keCorrHRTI, keSaturation, kdCorrYawDirection, MU, CTcorrected, CPcorrected, axialIndProvided, useWakeAngle, bd, useaUbU, aU, bU, adjustInitialWakeDiamToYaw, FLORISoriginal, cos_spread


# CALCULATING POWER COEFFICIENTS BASED ON POWER CURVE
def CpCt_solve(rotorDiameter, velocitiesTurbines, rpm, powercurve):
    Rtip = rotorDiameter/2.

    Cp_vec = np.zeros(np.size(rpm))
    axind_vec = np.zeros(np.size(rpm))
    tsr_rot = np.zeros(np.size(rpm))
    for i in range(np.size(rpm)):
        tsr_rot[i] = (rpm[i]*pi/30)*Rtip[i]/velocitiesTurbines[i]
        Cp_vec[i] = powercurve(tsr_rot[i])
        axind_vec[i] = fsolve(lambda a: (4.*a*(1-a)**2 - Cp_vec[i]), .3)

    Ct_vec = 4.*axind_vec*(1.-axind_vec)

    return Cp_vec, Ct_vec, axind_vec


# POWER CALCULATION FROM FLORIS MODEL
def floris_power(turbineXw, turbineYw, rotorDiameter, Vinf, rho, generator_efficiency, yaw_deg, rpm, powercurve):
    global power_corr

    pP, ke, keCorrArray, keCorrCT, Region2CT, kd, me, initialWakeDisplacement, initialWakeAngle, baselineCT, keCorrTI, baselineTI, keCorrHR, keCorrHRTI, keSaturation, kdCorrYawDirection, MU, CTcorrected, CPcorrected, axialIndProvided, useWakeAngle, bd, useaUbU, aU, bU, adjustInitialWakeDiamToYaw, FLORISoriginal, cos_spread = floris_parameters()

    p_near0 = 1.0

    Vin = np.ones_like(turbineXw)*Vinf

    # FOR REFERENCE:
    # axialInduction = np.ones_like(turbineXw)*1.0/3.0
    # Ct = 4.0*axialInduction*(1.0-axialInduction)
    # Cp = np.ones_like(turbineXw)*(0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2))

    Cp, Ct, axialInduction = CpCt_solve(rotorDiameter, Vin, rpm, powercurve)

    # Converging the velocities to match the power curve
    while True:
        velocitiesTurbines, _, _, _ = _florisunify.floris_unified(turbineXw, turbineYw, yaw_deg, rotorDiameter, Vinf, Ct, axialInduction, ke, kd, me, initialWakeDisplacement, bd, MU, aU, bU, initialWakeAngle, cos_spread, keCorrCT, Region2CT, keCorrArray, useWakeAngle, adjustInitialWakeDiamToYaw, axialIndProvided, useaUbU)

        if max(fabs(velocitiesTurbines-Vin)) <= 1e-8:
            break
        else:
            Vin = np.copy(velocitiesTurbines)
            Cp, Ct, axialInduction = CpCt_solve(rotorDiameter, Vin, rpm, powercurve)

    wt_power = generator_efficiency*(0.5*rho*((pi/4)*rotorDiameter**2)*Cp*velocitiesTurbines**3)/1000.
    power = sum(wt_power/power_corr)

    return velocitiesTurbines, wt_power/power_corr, power


# SPL CALCULATION BASED ON BPM ACOUSTIC MODEL
def bpm_noise(turbineX, turbineY, obs, rpm, blade, high, rad, c, alpha, velocitiesTurbines, windroseDirections):
    global noise_corr

    parallel_run = True # use joblib module for faster calculations using multiple processors
    # parallel_run = False

    nobs = np.size(obs[:,0])
    nwind = np.size(windroseDirections)

    c1 = c*0.25
    nu = 1.78e-5
    c0 = 343.2
    psi = 14.0
    AR = 17.

    if parallel_run == True:
        SPLf = Parallel(n_jobs=-1)(delayed(bpmnoise)(turbineX,turbineY,obs[j],windroseDirections[i],velocitiesTurbines[i],rpm[i],blade,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr) for i in range(nwind) for j in range(nobs) )
        SPL = np.array(SPLf)
    elif parallel_run == False:
        SPL = np.zeros(nobs*nwind) # setting up vector for the constraints
        k = 0
        for i in range(nwind):
            for j in range(nobs):
                SPL[k] = _bpmacoustic.turbinepos(turbineX,turbineY,obs[j],windroseDirections[i],velocitiesTurbines[i],rpm[i],blade,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr)
                k += 1

    return SPL


# SPL CALCULATION BASED ON BPM ACOUSTIC MODEL (used for joblib)
def bpmnoise(turbineX,turbineY,obs,windroseDirections,velocitiesTurbines,rpm,blade,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr):
    return _bpmacoustic.turbinepos(turbineX,turbineY,obs,windroseDirections,velocitiesTurbines,rpm,blade,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr)


# TURBINE SEPARATION FUNCTION
def sep_func(loc):
    global rotorDiameter

    space = 2. # rotor diameters apart

    n = np.size(loc)/2
    x = loc[0:n]
    y = loc[n:]
    sep = np.zeros((n-1)*n/2)

    k = 0
    for i in range(0, n):
        for j in range(i+1, n):
            sep[k] = (x[j]-x[i])**2+(y[j]-y[i])**2
            k += 1

    return sep - (space*rotorDiameter[0])**2


def obj_func(xdict):
    global rotorDiameter
    global Cp
    global Ct
    global axialInduction
    global generator_efficiency
    global yaw
    global rho
    global Vinf
    global windroseDirections
    global windFrequencies
    global powercurve
    global funcs
    global power_max
    global power_corr
    global obs
    global blade
    global high
    global rad
    global c
    global alpha

    x = xdict['xvars']*100. # rescaling back to original value
    y = xdict['yvars']*100. # rescaling back to original value
    rpm = xdict['rpm']*10. # rescaling back to original value
    funcs = {}

    nturb = np.size(x)
    nwind = np.size(windroseDirections)

    rpmw = np.zeros((nwind,nturb))
    k = 0
    for i in range(nwind):
        for j in range(nturb):
            rpmw[i,j] = rpm[k]
            k += 1

    # Calculate the average power production of the wind farm
    power_dir = np.zeros(nwind)
    veleff = np.zeros((nwind,nturb))
    windDirections = np.zeros_like(windroseDirections)
    for d in range(0, nwind):
        # Adjusting coordinate system
        windDirections[d] = 270. - windroseDirections[d]
        if windDirections[d] < 0.:
            windDirections[d] += 360.
        windDirectionRad = pi*windDirections[d]/180.0
        xw = x*cos(-windDirectionRad) - y*sin(-windDirectionRad)
        yw = x*sin(-windDirectionRad) + y*cos(-windDirectionRad)

        # Effective velocity in front of each turbine and power from each turbines
        veleff[d], _, power = floris_power(xw, yw, rotorDiameter, Vinf, rho, generator_efficiency, yaw[d], rpmw[d], powercurve)
        power_dir[d] = power*windFrequencies[d] # multiplying power by probability in that wind direction
    APP = sum(power_dir) # average power production

    funcs['obj'] = (-1.*APP)/1.0e4 #scaling to a unit of 1

    # Calculate the sound pressure level at each observer location
    SPL = bpm_noise(x, y, obs, rpmw, blade, high, rad, c, alpha, veleff, windroseDirections)
    funcs['SPL'] = (SPL)/10. #scaling to a unit of 1

    print 'Average power production:',APP,'kW ('+str(power_max*nturb)+' kW max)   Max SPL:',max(SPL),'dB'

    # Calculating the separation between each turbine
    sep = sep_func(np.append(x,y))
    funcs['sep'] = sep/100. #scaling to a unit of 1
    funcs['veleff'] = veleff/10. #scaling to a unit of 1

    fail = False

    return funcs, fail


# Main
if __name__ == "__main__":

    global rotorDiameter
    global Cp
    global Ct
    global axialInduction
    global generator_efficiency
    global power_max
    global power_corr
    global yaw
    global rho
    global Vinf
    global windroseDirections
    global windFrequencies
    global powercurve
    global funcs
    global rpm
    global obs
    global blade
    global high
    global rad
    global c
    global alpha
    global noise_corr

    SPLlim = float(argv[1]) # read in SPL limit
    # SPLlim = 100. # hardcode SPL limit

    # Wind farm data
    farmname = 'Lissett'
    # farmname = 'Rosiere'

    if farmname == 'Lissett':
        turbineX = np.array([460., 850., 320., 750., 180., 630., 100., 540., 110., 470., 100., 440.]) # x positions (m)
        turbineY = np.array([140., 320., 350., 500., 590., 770., 850., 1000., 1150., 1300., 1430., 1570.]) # y positions (m)
        obs = np.array([[1980.,-180.,2.],[2100.,1300.,2.],[1400.,2450.,2.],[80.,2660.,2.],[-2000.,1300.,2.],[-950.,1400.,2.],[-1250.,300.,2.]]) # observer positions (m)
        windFrequencies = np.array([0.0430,0.0362,0.0336,0.0247,0.0286,0.0122,0.0103,0.0407,0.1702,0.2265,0.1246,0.0909,0.0875,0.0279,0.0197,0.0234]) # percent of time the wind blows in each direction starting from North
        rho = 1.229 # air density (kg/m^3)
        Rhub = 0.9 # radial distance from the center to the rotor (m)
        high = 80. # height of turbine (m)
        rotor_diameter = 90. # diameter of turbine (m)
        rpm_max = 16.1 # maximum rotation rate (RPM)
        rpm_min = 0.#10.3 # minimum rotation rate (RPM)
        Vinf = 14. # free stream velocity (m/s)
        chord_corr = 1.645615 # correction factor for chord length
        power_corr = 1.842414141559 # correction factor for maximum power production
        power_max = 2500. # maximum power production of a single turbine (Nordex N90-2.5MW; kW)
        noise_corr = 0.8186289018703674 # correction factor for the SPL calculation

        # Wind farm boundaries
        xlow = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
        xupp = np.array([1200.,1200.,1200.,1200.,1200.,1200.,1200.,1200.,1200.,1200.,1200.,1200.])
        ylow = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
        yupp = np.array([1650.,1650.,1650.,1650.,1650.,1650.,1650.,1650.,1650.,1650.,1650.,1650.])

    if farmname == 'Rosiere':
        turbineX = np.array([1550., 1900., 2300., 2200., 2400., 2300., 2300., 2300., 220., 350., 350., 100., 900., 900., 800., 1100., 1300.]) # x positions (m)
        turbineY = np.array([4020., 4100., 3900., 3700., 3700., 3400., 3200., 3100., 1100., 1000., 800., 600., 1000., 900., 600., 300., 200.]) # y positions (m)
        obs = np.array([[1200.,4393.,2.],[1450.,4393.,2.],[2150.,4030.,2.],[2100.,4000.,2.],[1450.,3600.,2.],[1380.,3570.,2.],[600.,910.,2.],[0.,860.,2.],[0.,460.,2.],[500.,560.,2.],[900.,450.,2.],[1240.,440.,2.]]) # observer positions (m)
        windFrequencies = np.array([0.027185,0.044705,0.036813,0.035993,0.040072,0.053246,0.072519,0.096807,0.084873,0.094729,0.104432,0.126928,0.072634,0.040053,0.036927,0.032084]) # percent of time the wind blows in each direction starting from North
        rho = 1.176 # air density (kg/m^3)
        Rhub = 0.8382 # radial distance from the center to the rotor (m)
        high = 25. # height of turbine (m)
        rotor_diameter = 47. # diameter of turbine (m)
        rpm_max = 28.5 # maximum rotation rate (RPM)
        rpm_min = 0.#18. # minimum rotation rate (RPM)
        Vinf = 15. # free stream velocity (m/s)
        chord_corr = 2.190491 # correction factor for chord length
        power_corr = 2.08647545446 # correction factor for maximum power production
        power_max = 660. # maximum power production of a single turbine (Vestas V47-660kW; kW)
        noise_corr = 0.8695408271411205 # correction factor for the SPL calculation

        # Wind farm boundaries
        xlow = np.array([1404., 1734., 2121., 2121., 2121., 2121., 2121., 2121., 0., 0., 0., 0., 570., 570., 570., 1000., 1000.])
        xupp = np.array([1734., 2121., 2538., 2538., 2538., 2538., 2361., 2361., 570., 570., 570., 570., 1000., 1000., 1000., 1410., 1410.])
        ylow = np.array([4013., 4013., 3353., 3353., 3353., 3353., 3023., 3023., 400., 400., 400., 400., 810., 810., 400., 0., 0.])
        yupp = np.array([4393., 4393., 4013., 4013., 4013., 4013., 3353., 3353., 1410., 1410., 1410., 1410., 1210., 1210., 810., 400., 400.])

    nturb = np.size(turbineX)
    nwind = np.size(windFrequencies)
    nobs = np.size(obs[:,0])

    # Setting up wind rose
    increment = 360./nwind
    windroseDirections = np.zeros(nwind)
    for i in range(1, nwind):
        windroseDirections[i] = windroseDirections[i-1] + increment

    # Setting minimum and maximum rotation rates
    rpmlow = np.ones(nturb*nwind)*rpm_min
    rpmupp = np.ones(nturb*nwind)*rpm_max

    # Power Curve
    blade = 3. # number of turbine blades
    mu = 1.81206e-5 # fluid viscosity (kg/ms)

    # NREL 5 MW Turbine Geometry
    Rtip_nrel = 63.0 # turbine radius of NREL turbine
    Rtip = rotor_diameter/2. # turbine radius of respective wind farm
    r_nrel = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]) # radial positions (m)
    chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]) # chord lengths (m)
    alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]) # angle of attacks (deg)

    # Resize turbine geometry to sizes of turbine used on wind farm
    r_ratio = Rtip/Rtip_nrel
    rad = r_nrel*r_ratio
    c = chord_nrel*(r_ratio*chord_corr)

    # Reading in airfoil geometries
    import os
    afinit = CCAirfoil.initFromAerodynFile  # just for shorthand
    base = path.join(path.dirname(path.realpath('__file__')), 'CCBlade/test/5MW_AFFiles')
    basepath = base + os.path.sep

    # load all airfoils
    airfoil_types = [0]*8
    airfoil_types[0] = afinit(basepath + 'Cylinder1.dat')
    airfoil_types[1] = afinit(basepath + 'Cylinder2.dat')
    airfoil_types[2] = afinit(basepath + 'DU40_A17.dat')
    airfoil_types[3] = afinit(basepath + 'DU35_A17.dat')
    airfoil_types[4] = afinit(basepath + 'DU30_A17.dat')
    airfoil_types[5] = afinit(basepath + 'DU25_A17.dat')
    airfoil_types[6] = afinit(basepath + 'DU21_A17.dat')
    airfoil_types[7] = afinit(basepath + 'NACA64_A17.dat')

    # place at appropriate radial stations
    af_idx = [0, 0, 1, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7]

    af = [0]*len(rad)
    for i in range(len(rad)):
        af[i] = airfoil_types[af_idx[i]]

    # Aerodynamic properties of NREL turbine
    tilt = -5.0
    precone = 2.5
    yaw = 0.0
    shearExp = 0.2
    hubHt = high
    nSector = 8

    # create CCBlade object
    aeroanalysis = CCBlade(rad, c, alpha, af, Rhub, Rtip, blade, rho, mu, precone, tilt, yaw, shearExp, hubHt, nSector)

    tsr = np.linspace(0,20,100) # tip-speed ratio
    Uinf = Vinf*np.ones_like(tsr) # free stream wind speed
    Omega = ((Uinf*tsr)/Rtip)*(30./np.pi) # rotation rate (rad/s)
    pitch = np.ones_like(tsr)*0. # pitch (deg)

    # Calculating power coefficients at each tip-speed ratio
    CP,_,_ = aeroanalysis.evaluate(Uinf, Omega, pitch, coefficient=True)

    # Creating a power curve for the turbine (tip-speed ratio vs. power coefficient)
    powercurve = Akima1DInterpolator(tsr,CP)

    # Adjusting geometry for SPL calculations
    r_nrel = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333, 63.0]) # radial positions (m)
    rad = r_nrel*r_ratio

    # Initialize input variables
    rotorDiameter = np.ones(nturb)*rotor_diameter
    generator_efficiency = np.ones(nturb)*0.944
    yaw = np.ones((nwind,nturb))*0.
    rpm = np.ones(nwind*nturb)*rpm_max

    # Optimization
    optProb = Optimization('Wind_Farm_APP', obj_func)
    optProb.addObj('obj')

    # Design Variables (scaled to 1)
    nrpm = nturb*nwind
    optProb.addVarGroup('xvars', nturb, 'c', lower=xlow/100., upper=xupp/100., value=turbineX/100.) # x positions
    optProb.addVarGroup('yvars', nturb, 'c', lower=ylow/100., upper=yupp/100., value=turbineY/100.) # y positions
    optProb.addVarGroup('rpm', nrpm, 'c', lower=rpmlow/10., upper=rpmupp/10., value=rpm/10.) # rpm values

    # Constraints (scaled to 1)
    num_cons_sep = (nturb-1)*nturb/2
    optProb.addConGroup('sep', num_cons_sep, lower=0., upper=None) # separation between turbines
    num_cons_spl = nwind*nobs
    optProb.addConGroup('SPL', num_cons_spl, lower=0., upper=SPLlim/10.) # SPL limit

    opt = SNOPT()
    opt.setOption('Scale option',0)
    opt.setOption('Iterations limit',1000000)
    res = opt(optProb)
    # Printing optimization results (SNOPT format)
    print res

    # Final results (scaled back to original values)
    pow = np.array(-1*res.fStar)*1e4
    xf = res.xStar['xvars']*100.
    yf = res.xStar['yvars']*100.
    rpmf = res.xStar['rpm']*10.
    rpmfw = np.zeros((nwind,nturb))
    k = 0
    for i in range(nwind):
        for j in range(nturb):
            rpmfw[i,j] = rpmf[k]
            k += 1

    veleff = funcs['veleff']*10.
    SPL = funcs['SPL']*10.
    SPLw = np.zeros((nwind,nobs))
    k = 0
    for i in range(nwind):
        for j in range(nobs):
            SPLw[i,j] = SPL[k]
            k += 1

    # Printing optimization results (user-friendly format)
    print 'Wind Directions:',windroseDirections
    print 'APP:',pow,'kW'
    print 'X-locations (initial):',turbineX
    print 'X-locations (final):',xf
    print 'Y-locations (initial):',turbineY
    print 'Y-locations (final):',yf
    print 'RPM:',rpmfw # grouped by wind direction
    print 'Effective wind speeds:',veleff # grouped by wind direction
    print 'SPL:',SPLw # grouped by wind direction

    # Printing a text file containing the optimization results (saved into results folder)
    baseresults = path.join(path.dirname(path.realpath('__file__')), 'results')
    if farmname == 'Lissett':
        filename = baseresults + '/ParetoLissett/spl'+str(SPLlim)+'.txt'
    elif farmname == 'Rosiere':
        filename = baseresults + '/ParetoRosiere/spl'+str(SPLlim)+'.txt'
    target = open(filename,'w')
    if np.max(SPLw) > SPLlim:
        if np.max(SPLw) < (SPLlim+0.1):
            target.write('\nCONSTRAINT SATISFIED WITHIN LIMITS')
        else:
            target.write('\nCONSTRAINT EXCEEDED!')
    else:
        target.write('\nCONSTRAINT SATISFIED COMPLETELY')

    target.write('\nAverage Power Production: '+str(pow)+' kW\n')
    target.write('Max SPL: '+str(np.max(SPLw))+' dB\n')
    target.write('\nWind Directions: '+str(windroseDirections)+' degrees\n')
    target.write('X-locations (initial): '+str(turbineX)+' m\n')
    target.write('X-locations (final): '+str(xf)+' m\n')
    target.write('Y-locations (initial): '+str(turbineY)+' m\n')
    target.write('Y-locations (final): '+str(yf)+' m\n')
    target.write('RPM: '+str(rpmfw)+'\n')
    target.write('Effective wind speeds: '+str(veleff)+' m/s\n')
    target.write('SPL: '+str(SPLw)+' dB\n')
    target.close()
