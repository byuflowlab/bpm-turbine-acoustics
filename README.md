# BPM Turbine Acoustics

Turbine acoustic code using the BPM equations developed by Brooks, Pope, and Marcolini

Developed by Eric Tingey at Brigham Young University, 2015-2017

This code models the acoustic propagation of a wind turbine based on turbulent boundary layer edge noise, separation stall noise, tip vortex formation noise, laminar boundary layer vortex shedding noise, and trailing edge bluntness vortex shedding noise. Turbulent inflow noise is not assumed in this current code. The semi-empirical equations were developed from the NACA 0012 airfoil data and the blade segments used in the test file are based on the NREL 5-MW wind turbine. Scaling of the segments is based on the blade length specified.

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.


## Installation instructions

- system requirements: gfortran (using MinGW for Windows in order to use the commands here), python 2.7, numpy, scipy, matplotlib
- run:
```
python setup.py install
```
- or navigate to the directory and run the following command in the terminal to build the Fortran code:

Mac
```
$ cd wake_model
$ f2py -c  --opt=-O2 -m _bpmacoustic BPM_Acoustic_Model.f90
$ f2py -c  --opt=-O2 -m _bpmvawtacoustic BPM_VAWT_Acoustic_Model.f90
```

Windows
```
cd wake_model
python <\your\path\to\f2py.py> -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _bpmacoustic BPM_Acoustic_Model.f90
python <\your\path\to\f2py.py> -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _bpmvawtacoustic BPM_VAWT_Acoustic_Model.f90
```
(<\your\path\to\f2py.py>: most likely C:\Python27\Scripts\f2py.py)

## Running the Python code

This python code can be run from another file using:
```python
SPL_HAWT = _bpmacoustic.turbinepos(turbx, turby, obs, winddir, windvel, rpm, B, h, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr)
SPL_VAWT = _bpmvawtacoustic.turbinepos(ntheta, turbx, turby, obs, winddir, B, Hub, high, rad, c, c1, alpha, nu, c0, psi, AR, noise_corr, rot, velf, velx, vely, wakex, wakey)
```
