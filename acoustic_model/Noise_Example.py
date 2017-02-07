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

SPL_HAWT = _bpmacoustic.turbinepos(turbx, turby, obs, winddir, windvel, rpm, B, h, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)
Calculating the sound pressure level for a HAWT

Parameters
----------
turbx : array
    x-positions of all the turbines heard by an observer (east to west, meter)
turby : array
    y-positions of all the turbines heard by an observer (north to south, meter)
obs : array
    x-, y-, and z-position of a specified observer (E-W, N-S, height; meter)
winddir : float
    direction the wind blows from (180=N, 270=E, 0=S, 90=W; degree)
windvel : array
    wind velocity at each turbine in the specified wind direction (m/s)
rpm : array
    rotation rate of each turbine (RPM)
B : float
    number of blades on a turbine
h : float
    height of a turbine hub (meter)
rad : array
    radial positions of the blade geometry (meter)
c : array
    chord length at each radial segment (meter)
c1 : array
    distance from the pitch axis to leading edge at each radial segment (meter)
alpha : array
    angle of attack of each radial segment (degree)
nu : float
    kinematic viscosity of the air (m^2/s)
c0 : float
    speed of sound of the air (m/s)
psi : float
    solid angle of turbine blades between upper and lower sides of trailing edge (degree)
AR : float
    aspect ratio of turbine blades
noise_corr : float
    correction factor for SPL calculations (1=none, use if calculations differ from expected)

Returns
----------
SPL_HAWT : float
    sound pressure level calculated at observer location (dB)


SPL_VAWT = _bpmvawtacoustic.turbinepos(ntheta, turbx, turby, obs, winddir, B, Hub, high, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr, rot, velf, velx, vely, wakex, wakey)
Calculating the sound pressure level for a VAWT

Parameters
----------
ntheta : int
    number of points along blade flight path to calculate velocities
turbx : array
    x-positions of all the turbines heard by an observer (east to west, meter)
turby : array
    y-positions of all the turbines heard by an observer (north to south, meter)
obs : array
    x-, y-, and z-position of a specified observer (E-W, N-S, height; meter)
winddir : float
    direction the wind blows from (180=N, 270=E, 0=S, 90=W; degree)
B : float
    number of blades on a turbine
Hub : float
    hub height of a turbine (meter)
high : array
    height positions along the turbine blade (meter)
rad : float
    turbine radius (meter)
c : array
    chord length at each radial segment (meter)
c1 : array
    distance from the pitch axis to leading edge at each radial segment (meter)
alpha : array
    angle of attack of each radial segment (degree)
nu : float
    kinematic viscosity of the air (m^2/s)
c0 : float
    speed of sound of the air (m/s)
psi : float
    solid angle of turbine blades between upper and lower sides of trailing edge (degree)
AR : float
    aspect ratio of turbine blades
noise_corr : float
    correction factor for SPL calculations (1=none, use if calculations differ from expected)
rot : array
    rotation rate of each turbine (rad/s)
velf : float
    free stream wind speed (m/s)
velx : array
    the self-induced x-velocity of the turbine at each point along blade flight path (m/s)
vely : array
    the self-induced y-velocity of the turbine at each point along blade flight path (m/s)
wakex : array
    the wake influenced x-velcoity of the turbine at each point along the blade flight path (m/s)
wakey : array
    the wake influenced y-velcoity of the turbine at each point along the blade flight path (m/s)

Returns
----------
SPL_VAWT : float
    sound pressure level calculated at observer location (dB)
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
    if timeleft < 120.:
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
x_test = np.array([0.]) # x-locations of turbines (m)
y_test = np.array([0.]) # y-locations of turbines (m)
obs_test = np.array([0., 243.84, 0.]) # x-, y-, and z-location of the observer (m)
winddir_test = 0. # wind direction (deg)
rpm_test = np.array([28.5]) # rotation rate of the tubrines (rpm)
windvel_test = np.array([15.]) # wind velocity (m/s)
B_test = 3. # number of blades
h_test = 25. # height of the turbine hub (m)
noise_corr = 0.8697933840957954 # correction factor for noise

rad = np.array( [1.069324603174603, 2.088888888888889, 3.1084531746031745, 4.382936507936508, 5.912301587301587, 7.441666666666666, 8.971031746031747, 10.500396825396825, 12.029761904761905, 13.559126984126985, 15.088492063492065, 16.617857142857144, 18.147222222222222, 19.6765873015873, 20.951070634920637, 21.97063492063492, 22.990199206349207, 23.5] ) # radial positions (m)
c = np.array( [2.8941253867777776, 3.1490568155396828, 3.404805332214286, 3.7234696181666673, 3.8010929698730163, 3.6425779148095243, 3.4718065410555554, 3.2740712661825397, 3.062445496793651, 2.861441870269841, 2.660438243746032, 2.459434617222222, 2.258430990698413, 2.057427364174603, 1.889924342071429, 1.7044453858888888, 1.1594477481190477] ) # chord lengths (m)
alpha = np.array( [13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106] ) # angles of attack (deg)
c1 = c*0.25 # pitch axis (m)
AR = 17. # blade aspect ratio

nu = 1.78e-5 # kinematic viscosity (m^2/s)
c0 = 343.2 # speed of sound (m/s)
psi = 14.0 # solid angle (deg)

db_test_ros = _bpmacoustic.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, B_test, h_test, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)

print '\nTest Cases:'
print '\tRosiere Validation (47 dB): ', db_test_ros

##################################################################################
##################################################################################
##################################################################################

# SPL Test Paremeters (changed to whatever desired)
x_test = np.array([0.,5.]) # x-locations of turbines (m)
y_test = np.array([0.,5.]) # y-locations of turbines (m)
obs_test = np.array([0., 200., 0.]) # x-, y-, and z-location of the observer (m)
winddir_test = 0. # wind direction (deg)
rpm_test = np.array([28.5,28.5]) # rotation rate of the tubrines (rpm)
windvel_test = np.array([10.0,10.0]) # wind velocity (m/s)
B_test = 3. # number of blades
h_test = 25. # height of the turbine hub (m)
noise_corr = 0.8697933840957954 # correction factor for noise

rad = np.array( [1.069324603174603, 2.088888888888889, 3.1084531746031745, 4.382936507936508, 5.912301587301587, 7.441666666666666, 8.971031746031747, 10.500396825396825, 12.029761904761905, 13.559126984126985, 15.088492063492065, 16.617857142857144, 18.147222222222222, 19.6765873015873, 20.951070634920637, 21.97063492063492, 22.990199206349207, 23.5] ) # radial positions (m)
c = np.array( [2.8941253867777776, 3.1490568155396828, 3.404805332214286, 3.7234696181666673, 3.8010929698730163, 3.6425779148095243, 3.4718065410555554, 3.2740712661825397, 3.062445496793651, 2.861441870269841, 2.660438243746032, 2.459434617222222, 2.258430990698413, 2.057427364174603, 1.889924342071429, 1.7044453858888888, 1.1594477481190477] ) # chord lengths (m)
alpha = np.array( [13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106] ) # angles of attack (deg)
c1 = c*0.25 # pitch axis (m)
AR = 17. # blade aspect ratio

nu = 1.78e-5 # kinematic viscosity (m^2/s)
c0 = 343.2 # speed of sound (m/s)
psi = 14.0 # solid angle (deg)

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
        turbx = np.array([0.0]) # turbine x-positions (m)
        turby = np.array([0.0]) # turbine y-positions (m)
        winddir = 180.0 # wind direction (deg)
        rpm = np.array([16.1]) # rotation rates (RPM)
        windvel = np.array([13.5]) # wind velocities (m/s)
        B = 3. # number of blades
        h = 80. # hub height of turbine (m)

        rad = np.array( [0.8382, 2.047642857142857, 4.0, 5.952357142857142, 8.392857142857142, 11.321428571428571, 14.25, 17.17857142857143, 20.107142857142858, 23.035714285714285, 25.964285714285715, 28.892857142857146, 31.82142857142857, 34.75, 37.67857142857143, 40.11907142857143, 42.07142857142857, 44.023785714285715] ) # radial positions (m)
        c = np.array( [5.54194223, 6.030108795714287, 6.519839997857144, 7.130048205000001, 7.278688665714286, 6.97514919857143, 6.648140185, 6.269498169285715, 5.864257334285716, 5.479356772857143, 5.094456211428572, 4.70955565, 4.324655088571428, 3.939754527142857, 3.619004059285715, 3.26383159, 2.2202190921428575] ) # chord lengths (m)
        alpha = np.array( [13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106] ) # angles of attack (deg)
        c1 = c*0.25 # pitch axis (m)
        AR = 17. # blade aspect ratio

        nu = 1.78e-5 # kinematic viscosity (m^2/s)
        c0 = 343.2 # speed of sound (m/s)
        psi = 14.0 # solid angle (deg)

        # Setting up bounds of SPL calculations
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
        windcent = yp[0]/2.
        windlen = (yp[-1]-yp[0])/15.
        wx = windlen*sin(winddir*pi/180.)
        wy = windlen*cos(winddir*pi/180.)
        plt.annotate('', xy=(0-wx,windcent-wy), xycoords='data', xytext=(50.*sin(winddir*pi/180.),50.*cos(winddir*pi/180.)), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=20), fontsize=fs, color='k')
        plt.text(0., yp[0]/1.5, 'Wind', horizontalalignment='center', verticalalignment='top',fontsize=fs)

    elif turbine == 'vawt':
        ##################################################################################
        ##################################################################################
        ##################################################################################
        print '\nVertical-axis Wind Turbine Acoustic Plotting:'
        # Mariah Windspire 1.2 kW
        div = 5

        turbx = np.array([0.]) # turbine x-positions (m)
        turby = np.array([0.]) # turbine y-positions (m)
        velf = 8. # free stream wind speed (m/s)
        dia = 1.2 # turbine diameter (m)
        rad = dia/2. # turbine radius (m)
        tsrd = 2.625 # tip-speed ratio
        rot = np.ones_like(turbx)*tsrd*velf/rad # turbine rotation rates (rad/s)

        winddir = 180. # wind direction (deg)

        # Aerodynamic properties
        twist = 0.0
        delta = 0.0

        B = 3 # number of blades
        chord = 0.128 # chord length (m)
        c = np.ones(div)*chord # chord lengths over height positions (m)
        c1 = c*0.5 # pitch axis location (m)
        alpha = np.ones(div)*0.0 # angles of attack (deg)
        Hub = 2. # hub height (m)
        H = 6.1 # blade height (m)
        high = np.linspace(0,H,div+1) # height positions of the blade (m)
        AR = 5. # aspect ratio

        # Example self-induced wind velocities for isolated VAWT
        if rot[0] >= 0:
            velx = np.array([0.1300509758578449, 0.017234024387334474, -0.10055379878813524, -0.2061840120433015, -0.2858704822720968, -0.33119640487651264, -0.339089281964709, -0.3118508477141321, -0.25764478087546433, -0.1988684431703088, -0.17342520252077007, -0.15370972852456627, -0.1442486440082039, -0.15128511836734457, -0.13044595039768786, -0.08546908801023706, -0.02494466637533621, 0.03197120700568214, -0.011084310565649203, -0.16344583991429956, -0.3103017081750111, -0.41491602152175694, -0.4678517692849853, -0.4687091952147545, -0.4936006828598252, -0.531489515824587, -0.5750002772326674, -0.6644086256597977, -0.7427038690546394, -0.7791412301016253, -0.7634662797568859, -0.6910367446221855, -0.5625556663611492, -0.38570545877401596, -0.17512083613417861, 0.05394730308429136])
            vely = np.array([0.3070915729950489, 0.3483337599044662, 0.3456445783705026, 0.30465722308881105, 0.23435079208874923, 0.14686130947763262, 0.05584428582116426, -0.02466199807293788, -0.08071542678603759, -0.1040973357223225, -0.11184491158720485, -0.11999098742890327, -0.1290243531433885, -0.15323268231469964, -0.19116587377793248, -0.22288767900349615, -0.23658705088975251, -0.23010460906001576, -0.2122176454719008, -0.18721107170012113, -0.15682322144065844, -0.1274034585639046, -0.10121888330625133, -0.07291600358689668, -0.04069922131937081, -0.0084587366724128, 0.024310649639182044, 0.04943502865743335, 0.05726358249849845, 0.05178134955657134, 0.04222141788164915, 0.03878283394653364, 0.05048645736692045, 0.08399546577970662, 0.1430832510608585, 0.22711759171460444])
        elif rot[0] < 0:
            velx = np.array([0.03197122510517584, -0.024944642846871593, -0.08546913452518538, -0.13044596970446623, -0.15128513101776123, -0.14424860497342742, -0.1537097207976002, -0.17342525787701363, -0.19886843215397906, -0.2576446488126792, -0.31185084622227544, -0.33908921926319713, -0.33119639445512933, -0.2858705933082933, -0.20618403411859154, -0.10055381396353455, 0.017234011792231208, 0.13005100851350715, 0.05394731359530752, -0.17512087135369417, -0.38570545714893106, -0.5625557534839436, -0.6910369483427681, -0.763466224063827, -0.7791411750521885, -0.7427038148614863, -0.6644084421203552, -0.5750002331832405, -0.5314896507424663, -0.49360060975233755, -0.46870906862471234, -0.46785177276050205, -0.414916058862624, -0.31030181008786967, -0.16344575026009883, -0.011084268781767197])
            vely = np.array([0.23010459496417626, 0.2365870643937159, 0.22288768578535148, 0.19116585299835606, 0.1532326501818602, 0.12902433520071455, 0.11999099642974294, 0.11184488887236319, 0.10409726961483283, 0.08071541821126413, 0.024662027847672883, -0.05584425287347655, -0.14686123613799637, -0.23435076469016766, -0.3046572477295402, -0.34564459965523975, -0.3483337918857115, -0.30709160269194696, -0.2271175871718431, -0.14308322388898032, -0.08399546412807236, -0.05048645302266276, -0.03878281369863227, -0.042221435221775606, -0.05178134682427858, -0.05726357951893594, -0.04943502141218091, -0.024310615029862242, 0.00845877254564741, 0.040699256181776854, 0.07291600055261138, 0.10121886530845985, 0.12740344422184285, 0.15682322516209074, 0.18721106684450686, 0.2122176202642686])
        # Wake velocities from surronding turbines
        wakex = np.zeros_like(velx)
        wakey = np.zeros_like(vely)

        noise_corr = 1.0 # correction factor for noise

        # Setting up bounds of SPL calculations
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
        windcent = yp[0]/2.
        windlen = (yp[-1]-yp[0])/15.
        wx = windlen*sin(winddir*pi/180.)
        wy = windlen*cos(winddir*pi/180.)
        plt.annotate('', xy=(0-wx,windcent-wy), xycoords='data', xytext=(50.*sin(winddir*pi/180.),50.*cos(winddir*pi/180.)), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=20), fontsize=fs, color='k')
        plt.text(0., yp[0]/1.5, 'Wind', horizontalalignment='center', verticalalignment='top',fontsize=fs)

    plt.show()
