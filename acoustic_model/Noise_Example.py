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
import _bpmacoustic


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
noise_corr = -0.22243092339999748

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

db_test_ros = _bpmacoustic.turbinepos(x_test, y_test, obs_test, wind_test, rpm_test, windvel_test, B_test, h_test, rad, c, alpha, noise_corr)

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

db_test = _bpmacoustic.turbinepos(x_test, y_test, obs_test, wind_test, rpm_test, windvel_test, B_test, h_test, rad, c, alpha, noise_corr)

print 'Test SPL: ', db_test
