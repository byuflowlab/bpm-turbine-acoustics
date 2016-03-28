"""
BPM_complete.py

BPM Aeroacoustic Model by Eric Tingey, 2015
Based on acoustic theroy of Brooks, Pope & Marcolini

This code models the acoustic propagation of a wind turbine based on turbulent
boundary layer edge noise, separation stall noise, tip vortex formation noise,
laminar boundary layer vortex shedding noise, and trailing edge bluntness
vortex shedding noise. Turbulent inflow noise is not assumed in this current
code. The semi-empiracal equations were developed from the NACA 0012 airfoil
data and the blade segments used in the test file are based on the NREL 5-MW
wind turbine. Scaling of the segments is based on the blade length specified.
"""

import numpy as np
from math import factorial
from numpy import pi, log10, fabs, sqrt, sin, cos, arctan2, arctan
import matplotlib.pyplot as plt
import _bpmcomplete


def Dhfunc(theta_e, phi_e, M, Mc):
    """
    Directivity function for high-frequency noise (not high-angle separation or turbulent inflow noise); becomes inaccurate for theta_e approaching 180 deg
    """
    theta = theta_e * pi / 180
    phi = phi_e * pi / 180

    return (2 * (sin(theta / 2))**2 * (sin(phi))**2) / ((1 + M * cos(theta)) * (1 + (M - Mc) * cos(theta))**2)


def Dlfunc(theta_e, phi_e, M):
    """
    Directivity function for low-frequency noise
    """
    theta = theta_e * pi / 180
    phi = phi_e * pi / 180

    return ((sin(theta))**2 * (sin(phi))**2) / (1 + M * cos(theta))**4


def Afunc(ain, Re):
    """
    Spectral Function A
    Inputs:
    a: parameter = St/St_peak
    Re: Reynolds number of chord length
    """
    a = fabs(log10(ain))

    # Calculating Amin
    if(a < 0.204):
        Amin = sqrt(67.552 - 886.788 * a**2) - 8.22
    elif(a >= 0.204 and a <= 0.244):
        Amin = -32.665 * a + 3.981
    else:
        Amin = -142.79 * a**3 + 103.656 * a**2 - 57.757 * a + 6.006

    # Calculating Amax
    if(a < 0.13):
        Amax = sqrt(67.552 - 886.788 * a**2) - 8.219
    elif(a >= 0.13 and a <= 0.321):
        Amax = -15.901 * a + 1.098
    else:
        Amax = -4.669 * a**3 + 3.491 * a**2 - 16.699 * a + 1.149

    # Calculating a0
    if(Re < 9.52e4):
        a0 = 0.57
    elif(Re >= 9.52e4 and Re <= 8.57e5):
        a0 = -9.57e-13 * (Re - 8.57e5)**2 + 1.13
    else:
        a0 = 1.13

    # Calculating Amin(a0)
    if(a0 < 0.204):
        Amin0 = sqrt(67.552 - 886.788 * a0**2) - 8.22
    elif(a0 >= 0.204 and a0 <= 0.244):
        Amin0 = -32.665 * a0 + 3.981
    else:
        Amin0 = -142.79 * a0**3 + 103.656 * a0**2 - 57.757 * a0 + 6.006

    # Calculating Amax(a0)
    if(a0 < 0.13):
        Amax0 = sqrt(67.552 - 886.788 * a0**2) - 8.219
    elif(a0 >= 0.13 and a0 <= 0.321):
        Amax0 = -15.901 * a0 + 1.098
    else:
        Amax0 = -4.669 * a0**3 + 3.491 * a0**2 - 16.699 * a0 + 1.149

    AR = (-20.0 - Amin0) / (Amax0 - Amin0)

    return Amin + AR * (Amax - Amin)


def Bfunc(bin, Re):
    """
    Spectral Function B
    Inputs:
    b: parameter = Sts/St2
    Re: Reynolds number of chord length
    """
    b = fabs(log10(bin))

    # Calculating Bmin
    if(b < 0.13):
        Bmin = sqrt(16.888 - 886.788 * b**2) - 4.109
    elif(b >= 0.13 and b <= 0.145):
        Bmin = -83.607 * b + 8.138
    else:
        Bmin = -817.810 * b**3 + 355.210 * b**2 - 135.024 * b + 10.619

    # Calculating Bmax
    if(b < 0.10):
        Bmax = sqrt(16.888 - 886.788 * b**2) - 4.109
    elif(b >= 0.10 and b <= 0.187):
        Bmax = -31.330 * b + 1.854
    else:
        Bmax = -80.541 * b**3 + 44.174 * b**2 - 39.381 * b + 2.344

    # Calculating b0
    if(Re < 9.52e4):
        b0 = 0.30
    elif(Re >= 9.52e4 and Re <= 8.57e5):
        b0 = -4.48e-13 * (Re - 8.57e5)**2 + 0.56
    else:
        b0 = 0.56

    # Calculating Bmin(b0)
    if(b0 < 0.13):
        Bmin0 = sqrt(16.888 - 886.788 * b0**2) - 4.109
    elif(b0 >= 0.13 and b0 <= 0.145):
        Bmin0 = -83.607 * b0 + 8.138
    else:
        Bmin0 = -817.810 * b0**3 + 355.210 * b0**2 - 135.024 * b0 + 10.619

    # Calculating Bmax(b0)
    if(b0 < 0.10):
        Bmax0 = sqrt(16.888 - 886.788 * b**2) - 4.109
    elif(b0 >= 0.10 and b0 <= 0.187):
        Bmax0 = -31.330 * b + 1.854
    else:
        Bmax0 = -80.541 * b0**3 + 44.174 * b0**2 - 39.381 * b0 + 2.344

    BR = (-20 - Bmin0) / (Bmax0 - Bmin0)

    return Bmin + BR * (Bmax - Bmin)


def G1func(e):
    if(e <= 0.5974):
        G = 39.8 * log10(e) - 11.12
    elif(e <= 0.8545 and e > 0.5974):
        G = 98.409 * log10(e) + 2.0
    elif(e <= 1.17 and e > 0.8545):
        G = sqrt(2.484 - 506.25 * (log10(e))**2) - 5.076
    elif(e <= 1.674 and e > 1.17):
        G = -98.409 * log10(e) + 2.0
    else:
        G = -39.8 * log10(e) - 11.12

    return G


def G2func(d):
    if(d <= 0.3237):
        G = 77.852 * log10(d) + 15.328
    elif(d <= 0.5689 and d > 0.3237):
        G = 65.188 * log10(d) + 9.125
    elif(d <= 1.7579 and d > 0.5689):
        G = -114.052 * (log10(d))**2
    elif(d <= 3.0889 and d > 1.7579):
        G = -65.188 * log10(d) + 9.125
    else:
        G = -77.852 * log10(d) + 15.328

    return G


def G3func(alpha):

    return 171.04 - 3.03 * alpha


def G4func(hdav, psi):
    if(hdav <= 5.0):
        G = 17.5 * log10(hdav) + 157.5 - 1.114 * psi
    else:
        G = 169.7 - 1.114 * psi

    return G


def G5func(hdav, psi, StSt_peak):
    # finding G5 at phi = 14 deg
    eta = log10(StSt_peak)

    if(hdav < 0.25):
        mu = 0.1221
    elif(hdav < 0.62 and hdav >= 0.25):
        mu = -0.2175 * hdav + 0.1755
    elif(hdav < 1.15 and hdav >= 0.62):
        mu = -0.308 * hdav + 0.0596
    else:
        mu = 0.0242

    if(hdav <= 0.02):
        m = 0
    elif(hdav <= 0.5 and hdav > 0.02):
        m = 68.724 * (hdav) - 1.35
    elif(hdav <= 0.62 and hdav > 0.5):
        m = 308.475 * hdav - 121.23
    elif(hdav <= 1.15 and hdav > 0.62):
        m = 224.811 * hdav - 69.35
    elif(hdav <= 1.2 and hdav > 1.15):
        m = 1583.28 * hdav - 1631.59
    else:
        m = 268.344

    eta_0 = -sqrt((m**2 * mu**4) / (6.25 + m**2 * mu**2))
    k = 2.5 * sqrt(1 - (eta_0 / mu)**2) - 2.5 - m * eta_0

    if(eta < eta_0):
        G14 = m * eta + k
    elif(eta < 0 and eta >= eta_0):
        G14 = 2.5 * sqrt(1 - (eta / mu)**2) - 2.5
    elif(eta < 0.03616 and eta >= 0):
        G14 = sqrt(1.5625 - 1194.99 * eta**2) - 1.25
    else:
        G14 = -155.543 * eta + 4.375

    # finding G5 at psi = 0 deg
    hdav_prime = 6.724 * hdav**2 - 4.019 * hdav + 1.107

    eta0 = log10(StSt_peak)

    if(hdav_prime < 0.25):
        mu0 = 0.1221
    elif(hdav_prime < 0.62 and hdav_prime >= 0.25):
        mu0 = -0.2175 * hdav + 0.1755
    elif(hdav_prime < 1.15 and hdav_prime >= 0.62):
        mu0 = -0.308 * hdav + 0.0596
    else:
        mu0 = 0.0242

    if(hdav_prime <= 0.02):
        m0 = 0
    elif(hdav_prime <= 0.5 and hdav_prime > 0.02):
        m0 = 68.724 * (hdav_prime) - 1.35
    elif(hdav_prime <= 0.62 and hdav_prime > 0.5):
        m0 = 308.475 * hdav_prime - 121.23
    elif(hdav_prime <= 1.15 and hdav_prime > 0.62):
        m0 = 224.811 * hdav_prime - 69.35
    elif(hdav_prime <= 1.2 and hdav_prime > 1.15):
        m0 = 1583.28 * hdav_prime - 1631.59
    else:
        m0 = 268.344

    eta_00 = -sqrt((m0**2 * mu0**4) / (6.25 + m0**2 * mu0**2))
    k0 = 2.5 * sqrt(1 - (eta_00 / mu0)**2) - 2.5 - m0 * eta_00

    if(eta0 < eta_00):
        G0 = m0 * eta0 + k
    elif(eta0 < 0 and eta0 >= eta_00):
        G0 = 2.5 * sqrt(1 - (eta0 / mu0)**2) - 2.5
    elif(eta0 < 0.036 and eta0 >= 0):
        G0 = sqrt(1.5625 - 1194.99 * eta0**2) - 1.25
    else:
        G0 = -155.543 * eta0 + 4.375

    return G0 + 0.0714 * psi * (G14 - G0)


def TBLTE(f, V, L, c, r, theta_e, phi_e, alpha, nu, conv, trip):
    """
    Turbulent Boundary Layer Trailing Edge Noise
    Inputs:
    f: frequency of vortex shedding
    V: wind speed over the blade (m/s)
    L: span of the airfoil section (m)
    c: airfoil chord (m)
    r: effective observer distance (m)
    theta_e: angle of reference from the x-axis (position behind or in front of blade; 0-90 deg is behind)
    phi_e: angle of reference from the y-axis (position along the length of the blade; 0-90 deg is up)
    alpha: angle of attack (0 to 12.5 deg; greater is seperation stall)
    nu: kinematic viscosity (m^2/s)
    conv: convection factor for speed
    trip: boundary layer specification
    """
    c0 = 343.2  # speed of sound (m/s)
    M = V / c0
    Mc = conv * M
    Re = (V * c) / nu

    if trip == 'untrip':
        # UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
        d0 = c * (10**(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))**2))
        d0_d = c * (10**(3.0187 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))**2))
    elif trip == 'trip':
        # TRIPPED boundary layer at 0 deg- displacement thickness
        # this is from untripped, but since the pressure side doesn't matter,
        # use this value to carry on
        d0 = c * (10**(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))**2))
        if(Re <= 0.3e6):
            d0_d = c * 0.0601 * Re**(-0.114)
        else:
            d0_d = c * (10**(3.411 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))**2))

    # boundary layer on pressure side- thickness, displacement thickness
    dp = d0 * (10**(-0.04175 * alpha + 0.00106 * alpha**2))
    dp_d = d0_d * (10**(-0.0432 * alpha + 0.00113 * alpha**2))

    if trip == 'untrip':
        # UNTRIPPED boundary layer on suction side- displacement thickness
        if(alpha <= 7.5 and alpha >= 0):
            ds_d = d0_d * 10**(0.0679 * alpha)
        elif(alpha <= 12.5 and alpha > 7.5):
            ds_d = d0_d * 0.0162 * 10**(0.3066 * alpha)
        elif(alpha <= 25 and alpha > 12.5):
            ds_d = d0_d * 52.42 * 10**(0.0258 * alpha)
    elif trip == 'trip':
        # TRIPPED boundary layer on suction side- displacement thickness
        if(alpha <= 5 and alpha >= 0):
            ds_d = d0_d * 10**(0.0679 * alpha)
        elif(alpha <= 12.5 and alpha > 5):
            ds_d = d0_d * 0.381 * 10**(0.1516 * alpha)
        elif(alpha <= 25 and alpha > 12.5):
            ds_d = d0_d * 14.296 * 10**(0.0258 * alpha)

    Dh = Dhfunc(theta_e, phi_e, M, Mc)
    Dl = Dlfunc(theta_e, phi_e, M)

    Stp = (f * dp_d) / V
    Sts = (f * ds_d) / V

    St1 = 0.02 * M**(-0.6)

    if(alpha < 1.33):
        St2 = 1.0
    elif(alpha < 12.5 and alpha >= 1.33):
        St2 = 10**(0.0054 * (alpha - 1.33)**2)
    else:
        St2 = 4.72

    St_bar = 0.5 * (St1 + St2)

    apre = Stp / St1
    asuc = Sts / St_bar
    bang = Sts / St2

    gamma = 27.09 * M + 3.31
    gamma0 = 23.43 * M + 4.651
    beta = 72.65 * M + 10.74
    beta0 = -34.19 * M - 13.82

    if(Re < 2.5e5):
        K1 = -4.31 * log10(Re) + 156.3
    elif(Re >= 2.5e5 and Re <= 8.0e5):
        K1 = -9.0 * log10(Re) + 181.6
    else:
        K1 = 128.6

    if(alpha < (gamma0 - gamma)):
        K2 = K1 - 1000.0
    elif(alpha >= (gamma0 - gamma) and alpha <= (gamma0 + gamma)):
        K2 = K1 + sqrt(beta**2 - (beta / gamma)**2 * (alpha - gamma0)**2) + beta0
    else:
        K2 = K1 - 12.0

    Re_dp = (V * dp_d) / nu

    if(Re_dp <= 5000):
        DeltaK1 = alpha * (1.43 * log10(Re_dp)) - 5.29
    else:
        DeltaK1 = 0.0

    if(r < 1e-8):
        r = 1e-8

    if(alpha <= 12.5):
        Ap = Afunc(apre, Re)
        As = Afunc(asuc, Re)
        B = Bfunc(bang, Re)

        SPLp = 10 * log10((dp_d * M**5 * L * Dh) / r**2) + Ap + (K1 - 3.0) + DeltaK1
        SPLs = 10 * log10((dp_d * M**5 * L * Dh) / r**2) + As + (K1 - 3.0)
        SPLa = 10 * log10((ds_d * M**5 * L * Dh) / r**2) + B + K2

        return 10 * log10(10**(SPLp / 10) + 10**(SPLs / 10) + 10**(SPLa / 10))

    else:
        """
        Turbulent Boundary Layer Separation Stall Noise (TBLSS); this is were the airfoil is stalling and stall noise dominates
        SPLp = -infinity; 10**(SPLp/10) = 0
        SPLs = -infinity; 10**(SPLs/10) = 0
        """
        A = Afunc(bang, 3 * Re)

        SPLa = 10 * log10((ds_d * M**5 * L * Dl) / r**2) + A + K2

        return 10 * log10(10**(SPLa / 10))


def TBLTV(f, V, c, r, theta_e, phi_e, atip, conv, tip):
    """
    Turbulent Boundary Layer Tip Vortex Noise
    Inputs:
    f: frequency of vortex shedding
    V: wind speed over the blade (m/s)
    c: airfoil chord (m)
    r: effective observer distance (m)
    theta_e: angle of reference from the x-axis (position behind or in front of blade; 0-90 deg is behind)
    phi_e: angle of reference from the y-axis (position along the length of the blade; 0-90 deg is up)
    atip: angle of attack of the tip region (deg)
    conv: convection factor for speed
    tip: shape of the tip
    """
    c0 = 343.2 # speed of sound (m/s)
    M = V / c0
    Mc = conv * M

    Dh = Dhfunc(theta_e, phi_e, M, Mc)

    if tip == 'round':
        # rounded tip
        l = 0.008 * c * atip
    elif tip == 'flat':
        # flat tip
        if(atip <= 2 and atip >= 0):
            l = c * (0.0230 + 0.0169 * atip)
        else:
            l = c * (0.0378 + 0.0095 * atip)

    St = (f * l) / (V * (1 + 0.036 * atip))

    if(r < 1e-8):
        r = 1e-8
    return 10 * log10((M**5 * (1 + 0.036 * atip)**3 * l**2 * Dh) / r**2) - 30.5 * (log10(St) + 0.3)**2 + 126


def LBLVS(f, V, L, c, r, theta_e, phi_e, alpha, nu, conv):
    """
    Laminar Boundary Layer Vortex Shedding
    Inputs:
    f: frequency of vortex shedding
    V: wind speed (m/s)
    L: span of the airfoil section (m)
    c: airfoil chord (m)
    r: effective observer distance (m)
    theta_e: angle of reference from the x-axis (position behind or in front of blade; 0-90 deg is behind)
    phi_e: angle of reference from the y-axis (position along the length of the blade; 0-90 deg is up)
    alpha: angle of attack
    nu: kinematic viscosity (m^2/s)
    conv: convection factor for speed
    """
    c0 = 343.2 # speed of sound (m/s)
    M = V / c0
    Mc = conv * M

    Re = (V * c) / nu

    # UNTRIPPED boundary layer at 0 deg- thickness
    d0 = c * (10**(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))**2))
    # boundary layer on pressure side- thickness
    dp = d0 * (10**(-0.04175 * alpha + 0.00106 * alpha**2))

    St = (f * dp) / V

    Dh = Dhfunc(theta_e, phi_e, M, Mc)

    if(Re <= 1.3e5):
        St1 = 0.18
    elif(Re <= 4.0e5 and Re > 1.3e5):
        St1 = 0.001756 * Re**0.3931
    else:
        St1 = 0.28

    St_peak = St1 * 10**(-0.04 * alpha)

    e = St / St_peak

    G1 = G1func(e)

    if(alpha <= 3.0):
        Re0 = 10**(0.215 * alpha + 4.978)
    else:
        Re0 = 10**(0.12 * alpha + 5.263)

    d = Re / Re0

    G2 = G2func(d)
    G3 = G3func(alpha)

    if(r < 1e-8):
        r = 1e-8
    return 10 * log10((dp * M**5 * L * Dh) / r**2) + G1 + G2 + G3


def TEBVS(f, V, L, c, h, r, psi, theta_e, phi_e, alpha, nu, conv, trip):
    """
    Trailing Edge Bluntness Vortex Shedding Noise
    Inputs:
    f: frequency of vortex shedding
    V: wind speed over the blade (m/s)
    L: span of the airfoil section (m)
    c: airfoil chord (m)
    h: trailing edge thickness (m)
    r: effective observer distance (m)
    psi: solid angle between both airfoil surfaces just upstream of the trailing edge (deg)
    theta_e: angle of reference from the x-axis (position behind or in front of blade; 0-90 deg is behind)
    phi_e: angle of reference from the y-axis (position along the length of the blade; 0-90 deg is up)
    nu: kinematic viscosity (m^2/s)
    conv: convection factor for speed
    trip: boundary layer specification
    """
    c0 = 343.2 # speed of sound (m/s)
    M = V / c0
    Mc = conv * M

    Re = (V * c) / nu

    if trip == 'untrip':
        # UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
        d0 = c * (10**(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))**2))
        d0_d = c * (10**(3.0187 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))**2))
    elif trip == 'trip':
        # TRIPPED boundary layer at 0 deg- displacement thickness
        # this is from untripped, but since the pressure side doesn't matter,
        # use this value to carry on
        d0 = c * (10**(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))**2))
        if(Re <= 0.3e6):
            d0_d = c * 0.0601 * Re**(-0.114)
        else:
            d0_d = c * (10**(3.411 - 1.5397 * log10(Re) +
                             0.1059 * (log10(Re))**2))

    # boundary layer on pressure side- thickness, displacement thickness
    dp = d0 * (10**(-0.04175 * alpha + 0.00106 * alpha**2))
    dp_d = d0_d * (10**(-0.0432 * alpha + 0.00113 * alpha**2))

    if trip == 'untrip':
        # UNTRIPPED boundary layer on suction side- displacement thickness
        if(alpha <= 7.5 and alpha >= 0):
            ds_d = d0_d * 10**(0.0679 * alpha)
        elif(alpha <= 12.5 and alpha > 7.5):
            ds_d = d0_d * 0.0162 * 10**(0.3066 * alpha)
        elif(alpha <= 25 and alpha > 12.5):
            ds_d = d0_d * 52.42 * 10**(0.0258 * alpha)
    elif trip == 'trip':
        # TRIPPED boundary layer on suction side- displacement thickness
        if(alpha <= 5 and alpha >= 0):
            ds_d = d0_d * 10**(0.0679 * alpha)
        elif(alpha <= 12.5 and alpha > 5):
            ds_d = d0_d * 0.381 * 10**(0.1516 * alpha)
        elif(alpha <= 25 and alpha > 12.5):
            ds_d = d0_d * 14.296 * 10**(0.0258 * alpha)

    Dh = Dhfunc(theta_e, phi_e, M, Mc)
    St = (f * h) / V
    dav = (dp_d + ds_d) / 2

    hdav = h / dav

    if(hdav >= 0.2):
        St_peak = (0.212 - 0.0045 * psi) / (1 + 0.235 * (hdav)**-1 - 0.00132 * (hdav)**-2)
    else:
        St_peak = 0.1 * (hdav) + 0.095 - 0.00243 * psi

    StSt_peak = St / St_peak

    G4 = G4func(hdav, psi)
    G5 = G5func(hdav, psi, StSt_peak)

    if(r < 1e-8):
        r = 1e-8
    return 10 * log10((h * M**(5.5) * L * Dh) / r**2) + G4 + G5


def OASPL(r, theta_e, phi_e, rpm, Len, wind, B, Rhub, rad, c, alpha):
    """
    Computing the overall sound pressure level (OASPL) of a turbine defined below (in dB)
    Parameters:
    r: effective observer distance (m)
    theta_e: angle of reference from the x-axis (position behind or in front of blade; 0-90 deg is behind)
    phi_e: angle of reference from the y-axis (position along the length of the blade; 0-90 deg is up)
    rpm: rotational speed of the turbine (rpm)
    L: length of the turbine blade (m)
    wind: wind speed (m/s)
    B: number of turbine blades
    Rhub: the radial length of the hub section
    rad: radial positions along on the blade for segment division
    c: chord lengths at each radial section
    alpha: angles of attack at each radial section

    Returns:
    SPLoa: the overall sound pressure level at a given location from one turbine (dB)
    """
    # Using untripped or tripped boundary layer specficiation
    trip = 'untrip'
    # trip = 'trip'

    # Tip specfication
    tip = 'round'
    # tip = 'flat'

    nu = 1.78e-5  # kinematic viscosity of air (m^2/s)
    conv = 0.8  # convection factor for speed

    # Parameters of the wind turbine (f,V,L,c,h,alpha,atip)
    omega = (rpm * 2 * pi) / 60  # angular velocity (rad/sec)

    L = np.array([])
    V = np.array([])
    for i in range(np.size(rad)-1):
        wide = rad[i+1]-rad[i] #length of each radial section (m)
        rad_mid = (wide/2.)+rad[i] #midpoint of each radial section (m)
        vel = sqrt((omega*rad_mid)**2 + wind**2) #wind speed over the blade (m/s)
        L = np.append(L,wide)
        V = np.append(V,vel)

    h = 0.01 * c  # trailing edge thickness; 1% of chord length (m)
    atip = alpha[-1]  # angle of attack of the tip region (deg)
    psi = 1.520169109 #solid angle between both airfoil surfaces just upstream of the trailing edge (NACA 0012) (deg)

    # One-third octave band frequencies (Hz)
    f = np.array([10., 12.5, 16., 20., 25., 31.5, 40., 50., 63., 80., 100., 125., 160., 200., 250., 315., 400., 500., 630., 800., 1000., 1250., 1600., 2000., 2500., 3150., 4000., 5000., 6300., 8000., 10000., 12500., 16000., 20000.])

    nf = np.size(f)
    TE = np.linspace(0., 0., nf)
    TV = np.linspace(0., 0., nf)
    BLVS = np.linspace(0., 0., nf)
    BVS = np.linspace(0., 0., nf)
    SPLg = np.linspace(0., 0., nf)
    TE_temp = np.linspace(0., 0., np.size(L))
    TV_temp = np.linspace(0., 0., np.size(L))
    BLVS_temp = np.linspace(0., 0., np.size(L))
    BVS_temp = np.linspace(0., 0., np.size(L))
    refPres = 2e-5

    # Calculating sound pressure (Pa) for each noise source at each frequency and radial position
    for g in range(nf):
        for j in range(np.size(L)):
            TE_temp[j] = refPres * 10**(TBLTE(f[g],V[j],L[j],c[j],r,theta_e,phi_e,alpha[j],nu,conv,trip)/20)
            TV_temp[j] = refPres * 10**(TBLTV(f[g],V[j],c[j],r,theta_e,phi_e,atip,conv,tip)/20)
            BLVS_temp[j] = refPres * 10**(LBLVS(f[g],V[j],L[j],c[j],r,theta_e,phi_e,alpha[j],nu,conv)/20)
            BVS_temp[j] = refPres * 10**(TEBVS(f[g],V[j],L[j],c[j],h[j],r,psi,theta_e,phi_e,alpha[j],nu,conv,trip)/20)

        # Adding incoherent noises
        TE_t = np.sum(TE_temp)
        TV_t = np.sum(TV_temp)
        BLVS_t = np.sum(BLVS_temp)
        BVS_t = np.sum(BVS_temp)

        # Converting to sound pressure levels (SPL) (dB)
        TE[g] = 10 * log10(TE_t / (refPres)**2)
        TV[g] = 10 * log10(TV_t / (refPres)**2)
        BLVS[g] = 10 * log10(BLVS_t / (refPres)**2)
        BVS[g] = 10 * log10(BVS_t / (refPres)**2)
        SPLg[g] = 10 * log10(B*(10**(TE[g] / 10) + 10**(TV[g] / 10) + 10**(BLVS[g] / 10) + 10**(BVS[g] / 10))) # multiplied by number of blades (B)

    # A-weighting curve (dBA) for sound perception correction
    AdB = np.array([-70.4, -63.4, -56.7, -50.5, -44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1, -13.4, -10.9, -8.6, -6.6, -4.8, -3.2, -1.9, -0.8, 0.0, 0.6, 1.0, 1.2, 1.3, 1.2, 1.0, 0.5, -0.1, -1.1, -2.5, -4.3, -6.6, -9.3])
    SPLg = SPLg + AdB

    # Converting to sound pressure (Pa) for incoherent noise addition
    SPLp = refPres * 10**(SPLg / 20)

    # Converting back to SPL (dB)
    SPLoa = 10 * log10(np.sum(SPLp**2) / refPres**2)

    return SPLoa


def turbinepos(x, y, obs, wind, rpm, L, windvel, B, h, Rhub, rad, c, alpha, corr, py_for):
    """
    Placing a turbine in a specified location and finding the OASPL of the turbine with reference to an observer

    Parameters
    ----------
    x : array
        x locations of each of the turbines (x1,x2,...) (m)
    y : array
        y locations of each of the turbines (y1,y2,...) (m)
    obs : array
        the x, y, and z location of an observer (m)
    wind : float
        the direction of the wind flow in heading degrees (N=0,E=90,S=180,W=270)
    rpm : array
        rotational speed of each turbine (rpm)
    L : array
        length of the turbine blade (m)
    windvel : array
        wind speed at each turbine (m/s)
    B : float
        number of turbine blades
    h : float
        height of the turbine tower (m)
    Rhub : float
        the radial length of the hub section
    rad : array
        radial positions along on the blade for segment division (m)
    c : array
        chord lengths at each radial section (m)
    alpha : array
        angles of attack at each radial section (rad)
    corr : float
        SPL correction factor (dB)
    py_for : string
        specifying to use the Python ('python') or Fortran ('fortran') code

    Returns
    ----------
    turbinenoise : float
        the SPL at a specified observer location from all the turbines (dB)
    """
    if py_for == 'python':
        n = int(np.size(x))

        r = np.linspace(0, 0, n)
        theta_e = np.linspace(0, 0, n)
        phi_e = np.linspace(0, 0, n)
        t = np.linspace(0, 0, n)
        windrad = np.radians(wind)

        for i in range(n):
            r[i] = np.sqrt((obs[0] - x[i])**2 + (obs[1] - y[i])**2 + (obs[2] - h)**2)

            obsi = np.array([0.0, 0.0, 0.0])
            obsi[0] = obs[0] - x[i]
            obsi[1] = obs[1] - y[i]
            obsi[2] = obs[2]
            xi = 0.0
            yi = 0.0

            # Adjusting the coordinates for the wind direction
            rxy = np.sqrt((obsi[0])**2 + (obsi[1])**2)
            ang = arctan2(obsi[1], obsi[0]) + windrad

            obsi[0] = rxy * cos(ang)
            obsi[1] = rxy * sin(ang)

            phi_e[i] = fabs(np.degrees(arctan2((obsi[1] - yi), (obsi[0] - xi))))
            theta_e[i] = np.degrees(arctan2(fabs(h - obsi[2]), fabs(obsi[1] - yi)))
            if phi_e[i] < 1.0:
                phi_e[i] = 1.0
            if phi_e[i] > 179.0:
                phi_e[i] = 179.0
            # if phi_e[i] == 0.0 or phi_e[i] == 180.0:
            #     r[i] = r[i] - 0.2 #directivity adjustment based on work by Luis Vargas (Wind Turbine Noise Prediction)

        # Calculating the OASPL of each of the turbines in a for loop of "n" turbines (SPL function)
        for j in range(n):
            t[j] = OASPL(r[j], theta_e[j], phi_e[j], rpm[j], L, windvel[j], B, Rhub, rad, c, alpha)

        # Calculating the preliminary turbine_noise (10**(t1/10) + 10**(t2/10) + ...)
        turbine_noise = 10**(t[0] / 10)
        k = 1
        for num in range(n - 1):
            turbine_noise = turbine_noise + 10**(t[k] / 10)
            k = k + 1

        turbnoise = 10 * log10(turbine_noise)

        turbinenoise = turbnoise-corr #correction based on Rosiere valdiation (243.84 m, 800 ft should be 47 dB)

    elif py_for == 'fortran':
        turbinenoise = _bpmcomplete.turbinepos(x, y, obs, wind, rpm, L, windvel, B, h, rad, c, alpha, corr)

    return turbinenoise

# Running the Code
if __name__ == "__main__":

    # Rosiere Valdiation (243.84 m, 800 ft should be 47 dB;
    # http://www.mge.com/environment/green-power/wind/rosiere.htm)
    x_test = np.array([0.])
    y_test = np.array([0.])
    obs_test = np.array([0., 243.84, 0.])
    wind_test = 180.
    rpm_test = np.array([28.5])
    L_test = 23.5
    windvel_test = np.array([14.0])
    B_test = 3.
    h_test = 25.
    Rhub = 0.8382
    chord_corr = 2.190491

    # NREL 5-MW Turbine Specifications
    Rtip_nrel = 63.0

    r_nrel = np.array([1.5, 2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333])
    chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
    alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

    # Scaling the NREL turbine to the size of the Rosiere turbine
    r_ratio = L_test/Rtip_nrel
    rad = r_nrel*r_ratio
    rad[0] = Rhub
    c = chord_nrel*(r_ratio*chord_corr)

    db_test_ros = turbinepos(x_test, y_test, obs_test, wind_test, rpm_test, L_test, windvel_test, B_test, h_test, Rhub, rad, c, alpha, 30.781500225299993, 'fortran')

    print 'Rosiere Validation (47 dB): ', db_test_ros

    # SPL Test Paremeters (changed to whatever desired)
    x_test = np.array([0.,5.]) #x-location of turbine (m)
    y_test = np.array([0.,5.]) #y-location of turbine (m)
    obs_test = np.array([0., 200., 0.]) #x, y, and z-location of the observer (m)
    wind_test = 180. #wind direction (deg)
    rpm_test = np.array([28.5,28.5]) #roation rate of the tubrine (rpm)
    L_test = 23.4696 #length of the turbine blades (m)
    windvel_test = np.array([10.0,10.0]) #wind velocity (m/s)
    B_test = 3. #number of blades
    h_test = 25. #height of the turbine hub (m)
    Rhub = 0.8382 #radial length of the turbine hub (m)
    chord_corr = 2.190491

    # NREL 5-MW Turbine Specifications
    Rtip_nrel = 63.0

    r_nrel = np.array([1.5, 2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333])
    chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
    alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

    # Scaling the NREL turbine to the size of the specified turbine
    r_ratio = L_test/Rtip_nrel
    rad = r_nrel*r_ratio
    rad[0] = Rhub
    c = chord_nrel*(r_ratio*chord_corr)

    db_test = turbinepos(x_test, y_test, obs_test, wind_test, rpm_test, L_test, windvel_test, B_test, h_test, Rhub, rad, c, alpha, 30.781500225299993, 'fortran')

    print 'Test SPL: ', db_test
