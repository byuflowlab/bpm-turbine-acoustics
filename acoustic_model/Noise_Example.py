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
from numpy import sin,cos,pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import time,sys

import _bpmacoustic
import _bpmvawtacoustic

# Progress bar for plotting
def progress_bar(percent,tasks,runtime):
    bar_long = 40
    timeleft = (runtime)*(tasks*(1.-percent))
    if timeleft < 60.:
        status = 'Working... '+str(int(timeleft))+' seconds left'
    else:
        timeleft = timeleft/60.
        status = 'Working... '+str(int(timeleft))+' minutes left'
    if percent == 1:
        status = 'Complete\n'
    bar_seg = int(round(bar_long*percent))
    text = '\rStatus: [{0}] {1}% {2}'.format( '='*bar_seg + ' '*(bar_long - bar_seg), percent*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


# Option to plot a noise distribution around the turbines
plot_dist = True
# plot_dist = False # comment this out if desired on
turbine = 'hawt'
turbine = 'vawt'

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
windvel_test = np.array([15.])
B_test = 3.
h_test = 25.
Rhub = 0.8382
chord_corr = 2.190491
noise_corr = 0.8697933840957954

# NREL 5-MW Turbine Specifications
Rtip_nrel = 63.0

r_nrel = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333, 63.0])
chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

# Scaling the NREL turbine to the size of the Rosiere turbine
r_ratio = L_test/Rtip_nrel
rad = r_nrel*r_ratio
c = chord_nrel*(r_ratio*chord_corr)
c1 = c*0.25
AR = 17.

nu = 1.78e-5
c0 = 343.2
psi = 14.0

db_test_ros = _bpmacoustic.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, B_test, h_test, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)

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
noise_corr = 0.8697933840957954

nu = 1.78e-5
c0 = 343.2
psi = 14.0

# NREL 5-MW Turbine Specifications
Rtip_nrel = 63.0

r_nrel = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333, 63.0])
chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

# Scaling the NREL turbine to the size of the specified turbine
r_ratio = L_test/Rtip_nrel
rad = r_nrel*r_ratio
c = chord_nrel*(r_ratio*chord_corr)
c1 = c*0.25
AR = 17.

db_test = _bpmacoustic.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, B_test, h_test, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)

print '\tTest SPL: ', db_test

##################################################################################
##################################################################################
##################################################################################

if plot_dist == True:
    fs = 20 # font size
    fig = plt.figure(1,figsize=(8,7.8))
    fig.subplots_adjust(left=.12,right=.9,top=0.97,bottom=0.0)
    if turbine == 'hawt':
        ##################################################################################
        ##################################################################################
        ##################################################################################
        print '\nHorizontal-axis Wind Turbine Acoustic Plotting:'
        # Based on turbines found on Lissett Airfield Wind Farm
        # http://www.infinis.com/our-business/onshore-wind/our-operations/lissett-airfield/
        turbx = np.array([0.0])
        turby = np.array([0.0])
        winddir = 180.0
        rpm = np.array([16.1])
        windvel = np.array([13.5])
        B = 3.
        L = 45.
        h = 80.
        Rhub = 0.8382

        # NREL 5 MW Turbine
        Rtip_nrel = 63.0
        r_ratio = L/Rtip_nrel
        noise_corr = 0.8187768707500279

        r_nrel = np.array([1.5, 2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333])
        chord_nrel = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
        alpha = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

        # Scaling the NREL turbine to the size of the specified turbine
        rad = r_nrel*r_ratio
        rad[0] = Rhub
        chord_corr = 2.190491
        c = chord_nrel*r_ratio*chord_corr
        c1 = c*0.25
        AR = 17.
        nu = 1.78e-5
        c0 = 343.2
        psi = 14.0

        n = 100
        xp = np.linspace(-400.0, 400.0, n)
        yp = np.linspace(-400.0, 400.0, n)
        [X, Y] = np.meshgrid(xp, yp)
        F = np.zeros((n, n))

        point = 0
        time0 = time.time()
        for i in range(n):
            for j in range(n):
                F[i,j] = _bpmacoustic.turbinepos(turbx, turby, np.array([X[i,j],Y[i,j],0.]), winddir, windvel, rpm, B, h, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)
                point += 1
                runtime = time.time()-time0
                progress_bar(float(point)/(n*n),n*n,runtime)
                time0 = time.time()

        fig = plt.figure(1,figsize=(8,7.8))
        fig.subplots_adjust(left=.12,right=.9,top=0.97,bottom=0.0)

        lb = 20.0 #lower bound on velocity to display
        ub = 55.0 #upper bound on velocity to display
        ran = 100 #number of contours between the velocity bounds
        cbtix = np.linspace(lb,ub,8)
        bounds = np.linspace(lb,ub,ran)
        CS = plt.contourf(xp, yp, F, ran, cmap=plt.cm.viridis, vmax=ub, vmin=lb, levels=bounds)
        CB = plt.colorbar(CS,ticks=cbtix,orientation='horizontal',pad=0.13)
        CB.ax.tick_params(labelsize=fs)
        CB.ax.set_xlabel('Sound Pressure Level (dB)',fontsize=fs)
        ell1 = mpatches.Ellipse((L/2.,0.),L,L/5.625,color='k',fill=True)
        ell2 = mpatches.Ellipse((-L/2.,0.),L,L/5.625,color='k',fill=True)
        rect = mpatches.Rectangle((-L/11.25,-L/9.),L/5.625,L/3.,color='k',fill=True)
        plt.gca().add_patch(ell1)
        plt.gca().add_patch(ell2)
        plt.gca().add_patch(rect)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xlabel('Lateral Distance (m)',fontsize=fs)
        plt.ylabel('Downwind Distance (m)',fontsize=fs)
        plt.annotate('Wind', xy=(0,yp[0]/2.5),  xycoords='data', xytext=(0,-60), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=20), fontsize=fs,color='k')
    elif turbine == 'vawt':
        ##################################################################################
        ##################################################################################
        ##################################################################################
        print '\nVertical-axis Wind Turbine Acoustic Plotting:'
        # Mariah Windspire 1.2 kW
        div = 5

        turbx = np.array([0.])
        turby = np.array([0.])
        diat = np.ones_like(turbx)*1.2
        velf = 8.
        dia = 1.2
        rad = dia/2.
        tsrd = 2.625
        rot = np.ones_like(turbx)*tsrd*velf/rad

        winddir = 180.

        twist = 0.0
        delta = 0.0
        B = 3
        chord = 0.128
        c = np.ones(div)*chord
        c1 = c*0.5
        alpha = np.ones(div)*0.01
        Hub = 2.
        H = 6.1
        high = np.linspace(0,H,div+1)
        AR = 5.

        # Example self-induced wind velocities for isolated VAWT
        if rot >= 0:
            velx = np.array([0.1300509758578449, 0.017234024387334474, -0.10055379878813524, -0.2061840120433015, -0.2858704822720968, -0.33119640487651264, -0.339089281964709, -0.3118508477141321, -0.25764478087546433, -0.1988684431703088, -0.17342520252077007, -0.15370972852456627, -0.1442486440082039, -0.15128511836734457, -0.13044595039768786, -0.08546908801023706, -0.02494466637533621, 0.03197120700568214, -0.011084310565649203, -0.16344583991429956, -0.3103017081750111, -0.41491602152175694, -0.4678517692849853, -0.4687091952147545, -0.4936006828598252, -0.531489515824587, -0.5750002772326674, -0.6644086256597977, -0.7427038690546394, -0.7791412301016253, -0.7634662797568859, -0.6910367446221855, -0.5625556663611492, -0.38570545877401596, -0.17512083613417861, 0.05394730308429136])
            vely = np.array([0.3070915729950489, 0.3483337599044662, 0.3456445783705026, 0.30465722308881105, 0.23435079208874923, 0.14686130947763262, 0.05584428582116426, -0.02466199807293788, -0.08071542678603759, -0.1040973357223225, -0.11184491158720485, -0.11999098742890327, -0.1290243531433885, -0.15323268231469964, -0.19116587377793248, -0.22288767900349615, -0.23658705088975251, -0.23010460906001576, -0.2122176454719008, -0.18721107170012113, -0.15682322144065844, -0.1274034585639046, -0.10121888330625133, -0.07291600358689668, -0.04069922131937081, -0.0084587366724128, 0.024310649639182044, 0.04943502865743335, 0.05726358249849845, 0.05178134955657134, 0.04222141788164915, 0.03878283394653364, 0.05048645736692045, 0.08399546577970662, 0.1430832510608585, 0.22711759171460444])
        elif rot < 0:
            velx = np.array([0.03197122510517584, -0.024944642846871593, -0.08546913452518538, -0.13044596970446623, -0.15128513101776123, -0.14424860497342742, -0.1537097207976002, -0.17342525787701363, -0.19886843215397906, -0.2576446488126792, -0.31185084622227544, -0.33908921926319713, -0.33119639445512933, -0.2858705933082933, -0.20618403411859154, -0.10055381396353455, 0.017234011792231208, 0.13005100851350715, 0.05394731359530752, -0.17512087135369417, -0.38570545714893106, -0.5625557534839436, -0.6910369483427681, -0.763466224063827, -0.7791411750521885, -0.7427038148614863, -0.6644084421203552, -0.5750002331832405, -0.5314896507424663, -0.49360060975233755, -0.46870906862471234, -0.46785177276050205, -0.414916058862624, -0.31030181008786967, -0.16344575026009883, -0.011084268781767197])
            vely = np.array([0.23010459496417626, 0.2365870643937159, 0.22288768578535148, 0.19116585299835606, 0.1532326501818602, 0.12902433520071455, 0.11999099642974294, 0.11184488887236319, 0.10409726961483283, 0.08071541821126413, 0.024662027847672883, -0.05584425287347655, -0.14686123613799637, -0.23435076469016766, -0.3046572477295402, -0.34564459965523975, -0.3483337918857115, -0.30709160269194696, -0.2271175871718431, -0.14308322388898032, -0.08399546412807236, -0.05048645302266276, -0.03878281369863227, -0.042221435221775606, -0.05178134682427858, -0.05726357951893594, -0.04943502141218091, -0.024310615029862242, 0.00845877254564741, 0.040699256181776854, 0.07291600055261138, 0.10121886530845985, 0.12740344422184285, 0.15682322516209074, 0.18721106684450686, 0.2122176202642686])
        # Wake velocities from surronding turbines
        wakex = np.zeros_like(velx)
        wakey = np.zeros_like(vely)

        noise_corr = 1.0

        n = 100
        xp = np.linspace(-10.0, 10.0, n)
        yp = np.linspace(-10.0, 10.0, n)
        [X, Y] = np.meshgrid(xp, yp)
        F = np.zeros((n, n))

        point = 0
        time0 = time.time()
        for i in range(n):
            for j in range(n):
                time0 = time.time()
                F[i,j] = _bpmvawtacoustic.turbinepos(36, turbx, turby, np.array([X[i,j],Y[i,j],0.]), winddir, B, Hub, high, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr, rot, velf, velx, vely, wakex, wakey)
                point += 1
                runtime = time.time()-time0
                progress_bar(float(point)/(n*n),n*n,runtime)
                time0 = time.time()

        lb = 66.0 #lower bound on velocity to display
        ub = 75.0 #upper bound on velocity to display
        ran = 100 #number of contours between the velocity bounds
        cbtix = np.linspace(lb,ub,10)
        bounds = np.linspace(lb,ub,ran)
        CS = plt.contourf(xp, yp, F, ran, cmap=plt.cm.viridis, vmax=ub, vmin=lb, levels=bounds)
        CB = plt.colorbar(CS,ticks=cbtix,orientation='horizontal',pad=0.13)
        CB.ax.tick_params(labelsize=fs)
        CB.ax.set_xlabel('Sound Pressure Level (dB)',fontsize=fs)
        circ = mpatches.Circle((0.,0.),rad,color='k',linestyle=':',fill=False,linewidth=1)
        ell1 = mpatches.Ellipse((0.,rad),c[0],c[0]*0.25,color='k',fill=True)
        ell2 = mpatches.Ellipse((rad*cos(210.*pi/180.),rad*sin(210.*pi/180.)),c[0],c[0]*0.25,angle=120.,color='k',fill=True)
        ell3 = mpatches.Ellipse((rad*cos(330.*pi/180.),rad*sin(330.*pi/180.)),c[0],c[0]*0.25,angle=240.,color='k',fill=True)
        plt.gca().add_patch(circ)
        plt.gca().add_patch(ell1)
        plt.gca().add_patch(ell2)
        plt.gca().add_patch(ell3)
        plt.plot((0.,0.),(0.,rad),'k',linewidth=1.)
        plt.plot((0.,rad*cos(210.*pi/180.)),(0.,rad*sin(210.*pi/180.)),'k',linewidth=1.)
        plt.plot((0.,rad*cos(330.*pi/180.)),(0.,rad*sin(330.*pi/180.)),'k',linewidth=1.)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xlabel('Lateral Distance (m)',fontsize=fs)
        plt.ylabel('Downwind Distance (m)',fontsize=fs)
        plt.annotate('Wind', xy=(0,yp[0]/2.5),  xycoords='data', xytext=(0,-60), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=20), fontsize=fs,color='k')

    plt.show()
