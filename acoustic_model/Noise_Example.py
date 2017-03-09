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
        rot = -np.ones_like(turbx)*tsrd*velf/rad # turbine rotation rates (rad/s)

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
            velx = np.array([0.1581151375458653, 0.10817861642691516, 0.05271722269052151, -0.005636963356682027, -0.06462855401801865, -0.12211840606397469, -0.17606531327739475, -0.22454836072756534, -0.26591094099284257, -0.29892029852828994, -0.32271131518655893, -0.3367743339012266, -0.340977475870952, -0.33559057369676387, -0.32134037162312085, -0.29946352130742177, -0.27152929206141485, -0.23604607660586963, -0.2085377555432528, -0.19265927242501976, -0.1800813387705058, -0.16907889094051037, -0.1592405289043837, -0.15093195328836703, -0.14579999636734278, -0.15013200453879824, -0.1539065803686135, -0.14937706078292562, -0.13854522618171447, -0.1212428677070831, -0.09806261277625498, -0.07022738579625547, -0.039544956432850505, -0.008667045305260223, 0.019379834772031656, 0.04233982065419682, 0.02115782441374853, -0.04595328640008716, -0.12131367018123537, -0.19932953858582095, -0.27338072719294815, -0.3384928776200799, -0.3919205073377353, -0.4322890747188202, -0.4594284266942592, -0.4744193768966408, -0.474424411105666, -0.472690944301474, -0.48560970445396995, -0.5028869738676074, -0.5217852222233792, -0.5416987495580154, -0.5632422209425544, -0.5889690485790483, -0.6310420173932223, -0.683601555846162, -0.7240942781081103, -0.7548518030713258, -0.7739100945339932, -0.7797180013243837, -0.7712779452693823, -0.7480398196448479, -0.7098430378994942, -0.6569450155943053, -0.5900610363804705, -0.5104400701382934, -0.4199484267660619, -0.3207852043509741, -0.21516635388249686, -0.1051645475568925, 0.007492062133359142, 0.1213741580996115])
            vely = np.array([0.2954622495968541, 0.3303295632702701, 0.35286995770127133, 0.36391470352397404, 0.36386324524958497, 0.3531972623290628, 0.3326396785337387, 0.3032386979862731, 0.266392703336669, 0.22375475172703563, 0.1771260209591327, 0.12841649041760342, 0.07961038921870778, 0.03272558795360017, -0.010247549071908469, -0.04747505366580312, -0.07750434019245106, -0.09584067087189048, -0.10238848096730101, -0.10560392867173869, -0.10978250575763257, -0.11422551833186928, -0.11831277788726467, -0.12154621832900053, -0.12383857648811261, -0.13007105126052437, -0.14479047330868855, -0.16460349133930718, -0.18547147659398028, -0.2054998121868598, -0.22269754957066956, -0.23534919446728167, -0.2421446519347792, -0.2425598318389634, -0.23746862546204278, -0.2296548361184321, -0.22098818661245975, -0.2104650239183888, -0.1972989094934629, -0.18209498382609607, -0.16592479353613943, -0.14992151430321535, -0.13485909722550993, -0.12102092939733074, -0.10820761292888526, -0.09586799553996818, -0.08248621580451243, -0.0666065166364508, -0.049593080355833756, -0.03307694345314293, -0.016939302429485895, -0.0009224440219106195, 0.015327382622703712, 0.03226108080877977, 0.04772867877295569, 0.05710712207862905, 0.060198950306511255, 0.05930900235084881, 0.05524656719113472, 0.049235427342751055, 0.04271207375957689, 0.037155268972744154, 0.033965809742782835, 0.03438780808596107, 0.039474058658145975, 0.05008854210383916, 0.06689135331720122, 0.09031568108960523, 0.12056403933715584, 0.15757089328302187, 0.20087104731696714, 0.24939807485416365])
        elif rot[0] < 0:
            velx = np.array([0.042339820654128994, 0.01937983477196117, -0.008667045305325024, -0.03954495643292008, -0.07022738579632888, -0.09806261277631743, -0.12124286770699459, -0.13854522618176596, -0.14937706078299773, -0.15390658036815474, -0.15013200453884734, -0.14579999636729796, -0.1509319532883437, -0.15924052890440524, -0.16907889094054745, -0.18008133877056032, -0.19265927242516442, -0.20853775554259424, -0.23604607660451193, -0.27152929206162774, -0.2994635213072936, -0.3213403716230264, -0.3355905736967202, -0.3409774758709287, -0.33677433390123224, -0.32271131518656715, -0.2989202985283245, -0.26591094099288765, -0.22454836072762743, -0.17606531327745523, -0.12211840606405298, -0.06462855401810438, -0.0056369633567841014, 0.05271722269041364, 0.108178616426798, 0.15811513754571851, 0.1213741580994426, 0.007492062133239851, -0.10516454755700122, -0.2151663538825826, -0.32078520435105135, -0.4199484267661296, -0.5104400701383132, -0.5900610363805063, -0.6569450155943202, -0.7098430378994866, -0.7480398196448145, -0.7712779452693133, -0.7797180013242995, -0.7739100945338512, -0.754851803071098, -0.7240942781076514, -0.6836015558461201, -0.6310420173908146, -0.5889690485787269, -0.5632422209428043, -0.541698749558044, -0.5217852222233634, -0.5028869738675436, -0.4856097044538084, -0.47269094430123476, -0.4744244111054857, -0.474419376896299, -0.4594284266945807, -0.4322890747189185, -0.39192050733749667, -0.33849287762013286, -0.27338072719302153, -0.19932953858588537, -0.12131367018129491, -0.045953286400153115, 0.021157824413696527])
            vely = np.array([0.22965483611839996, 0.2374686254619812, 0.24255983183891067, 0.24214465193470305, 0.2353491944671765, 0.22269754957053006, 0.2054998121867485, 0.18547147659387145, 0.16460349133907587, 0.1447904733086219, 0.13007105126057877, 0.1238385764880771, 0.12154621832897021, 0.11831277788721038, 0.11422551833176055, 0.10978250575745536, 0.1056039286712956, 0.102388480966803, 0.09584067087215903, 0.07750434019297735, 0.04747505366599366, 0.0102475490721258, -0.032725587953395446, -0.07961038921851649, -0.1284164904174351, -0.17712602095897453, -0.22375475172688866, -0.2663927033365169, -0.3032386979861349, -0.332639678533619, -0.3531972623289444, -0.3638632452494882, -0.3639147035238786, -0.3528699577011961, -0.3303295632702049, -0.2954622495967796, -0.24939807485414872, -0.20087104731695632, -0.15757089328303264, -0.12056403933716436, -0.09031568108963323, -0.06689135331724125, -0.0500885421038725, -0.03947405865819211, -0.03438780808601541, -0.03396580974283986, -0.037155268972800096, -0.04271207375963295, -0.04923542734280532, -0.055246567191180394, -0.059309002350854595, -0.06019895030655426, -0.05710712207898175, -0.0477286787734903, -0.032261080808760646, -0.015327382622355489, 0.0009224440220844055, 0.016939302429559523, 0.03307694345317264, 0.04959308035581677, 0.06660651663639519, 0.08248621580420334, 0.09586799553978073, 0.10820761292912937, 0.12102092939751337, 0.1348590972256211, 0.1499215143032721, 0.16592479353618175, 0.18209498382611064, 0.19729890949347248, 0.21046502391836935, 0.22098818661240768])

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
                F[i,j] = _bpmvawtacoustic.turbinepos(np.size(velx), turbx, turby, np.array([X[i,j],Y[i,j],0.]), winddir, B, Hub, high, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr, rot, velf, velx, vely, wakex, wakey)
                point += 1
                runtime = time.time()-time0
                progress_bar(float(point)/(n*n),n*n,runtime)
                time0 = time.time()

        lb = 64.0 #lower bound on velocity to display
        ub = 73.0 #upper bound on velocity to display
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
