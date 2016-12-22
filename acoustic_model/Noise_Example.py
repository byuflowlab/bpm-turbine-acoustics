"""
Noise_Example.py

BPM Aeroacoustic Model by Eric Tingey, created in 2015
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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

import _bpmacoustic

# Option to plot a noise distribution around the turbines
plot_dist = True
# plot_dist = False # comment this out if desired on

##################################################################################
##################################################################################
##################################################################################

# Rosiere Valdiation (243.84 m, 800 ft should be 47 dB)
# http://www.mge.com/environment/green-power/wind/rosiere.htm
x_test = np.array([0.])
y_test = np.array([0.])
obs_test = np.array([0., 243.84, 0.])
winddir_test = 180.
rpm_test = np.array([28.5])
L_test = 23.5
windvel_test = np.array([14.0])
B_test = 3.
h_test = 25.
Rhub = 0.8382
chord_corr = 2.190491
noise_corr = 1.7862546970999986

# NREL 5-MW Turbine Specifications
Rtip_nrel = 63.0

r_nrel = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333, 63.0])
chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

# Scaling the NREL turbine to the size of the Rosiere turbine
r_ratio = L_test/Rtip_nrel
rad = r_nrel*r_ratio
c = chord_nrel*(r_ratio*chord_corr)

db_test_ros = _bpmacoustic.turbinepos(x_test, y_test, obs_test, winddir_test, rpm_test, windvel_test, B_test, h_test, rad, c, alpha, noise_corr)

if plot_dist == False:
    print '\nTest Cases:'
    print '\tRosiere Validation (47 dB): ', db_test_ros

##################################################################################
##################################################################################
##################################################################################

# SPL Test Paremeters (changed to whatever desired)
x_test = np.array([0.,5.]) #x-locations of turbines (m)
y_test = np.array([0.,5.]) #y-locations of turbines (m)
obs_test = np.array([0., 200., 0.]) #x, y, and z-location of the observer (m)
winddir_test = 180. #wind direction (deg)
rpm_test = np.array([28.5,28.5]) #roation rate of the tubrine (rpm)
L_test = 23.4696 #length of the turbine blades (m)
windvel_test = np.array([10.0,10.0]) #wind velocity (m/s)
B_test = 3. #number of blades
h_test = 25. #height of the turbine hub (m)
Rhub = 0.8382 #radial length of the turbine hub (m)
chord_corr = 2.190491

# NREL 5-MW Turbine Specifications
Rtip_nrel = 63.0

r_nrel = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333, 63.0])
chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

# Scaling the NREL turbine to the size of the specified turbine
r_ratio = L_test/Rtip_nrel
rad = r_nrel*r_ratio
c = chord_nrel*(r_ratio*chord_corr)

db_test = _bpmacoustic.turbinepos(x_test, y_test, obs_test, winddir_test, rpm_test, windvel_test, B_test, h_test, rad, c, alpha, noise_corr)

if plot_dist == False:
    print '\tTest SPL: ', db_test

##################################################################################
##################################################################################
##################################################################################

if plot_dist == True:
    # Based on turbines found on Lissett Airfield Wind Farm
    # http://www.infinis.com/our-business/onshore-wind/our-operations/lissett-airfield/
    turbx = np.array([0.0])
    turby = np.array([0.0])
    winddir = 0.0
    rpm = np.array([16.1])
    windvel = np.array([13.5])
    B = 3.
    L = 45.
    h = 80.
    Rhub = 0.8382

    # NREL 5 MW Turbine
    Rtip_nrel = 63.0
    r_ratio = L/Rtip_nrel
    noise_corr = 1.7862546970999986 #noise correction to match the Rosiere Wind Farm SPL of 47 dB

    r_nrel = np.array([1.5, 2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333])
    chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
    alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

    # Scaling the NREL turbine to the size of the specified turbine
    rad = r_nrel*r_ratio
    rad[0] = Rhub
    chord_corr = 2.190491
    c = chord_nrel*r_ratio*chord_corr

    n = 50
    xp = np.linspace(-250.0, 250.0, n)
    yp = np.linspace(-250.0, 250.0, n)
    [X, Y] = np.meshgrid(xp, yp)
    F = np.zeros((n, n))

    point = 0
    for i in range(n):
        for j in range(n):
            F[i,j] = _bpmacoustic.turbinepos(turbx, turby, np.array([X[i,j],Y[i,j],0.]), winddir, rpm, windvel, B, h, rad, c, alpha, noise_corr)
            point += 1
            print 'Calculating', point, 'of', n**2

    fig = plt.figure(1,figsize=(8,7.8))
    fig.subplots_adjust(left=.12,right=.9,top=0.97,bottom=0.0)

    lb = 46.0 #lower bound on velocity to display
    ub = 86.0 #upper bound on velocity to display
    ran = 200 #number of contours between the velocity bounds
    cbtix = np.linspace(lb,ub,6)
    bounds = np.linspace(lb,ub,ran)
    CS = plt.contourf(xp, yp, F, ran, cmap=plt.cm.parula, vmax=ub, vmin=lb, levels=bounds)
    # CSn = plt.contour(xp, yp, F, cmap=plt.cm.parula)
    CB = plt.colorbar(CS,ticks=cbtix,orientation='horizontal',pad=0.13)
    CB.ax.tick_params(labelsize=15)
    CB.ax.set_xlabel('Sound Pressure Level (dB)',fontsize=16)
    ell1 = mpatches.Ellipse((L/2.,0.),L,8.,color='k',fill=True)
    ell2 = mpatches.Ellipse((-L/2.,0.),L,8.,color='k',fill=True)
    rect = mpatches.Rectangle((-4.,-5.),8,15,color='k',fill=True)
    plt.gca().add_patch(ell1)
    plt.gca().add_patch(ell2)
    plt.gca().add_patch(rect)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Lateral Distance (m)',fontsize=16)
    plt.ylabel('Downwind Distance (m)',fontsize=16)
    plt.annotate('Wind', xy=(0,-130),  xycoords='data', xytext=(0,-60), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=20), fontsize=14,color='k')

    print '\nTest Cases:'
    print '\tRosiere Validation (47 dB): ', db_test_ros
    print '\tTest SPL: ', db_test

    plt.show()
