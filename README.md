# BPM Turbine Acoustics

Turbine acoustic code using the BPM equations developed by Brooks, Pope, and Marcolini

Developed by Eric Tingey at Brigham Young University, 2015

This code models the acoustic propagation of a wind turbine based on turbulent boundary layer edge noise, separation stall noise, tip vortex formation noise, laminar boundary layer vortex shedding noise, and trailing edge bluntness vortex shedding noise. Turbulent inflow noise is not assumed in this current code. The semi-empiracal equations were developed from the NACA 0012 airfoil data and the blade segments used in the test file are based on the NREL 5-MW wind turbine. Scaling of the segments is based on the blade length specified.

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.


## Installation instructions

- system requirements: gfortran (using MinGW for Windows in order to use the commands here), python 2.7, numpy, scipy
- navigate to the directory and run the following command in the terminal to build the Fortran code:

Mac
```
$ cd wake_model
$ f2py -c  --opt=-O2 -m _bpmcomplete BPM_complete.f90
```

Windows
```
cd wake_model
python <\your\path\to\f2py.py> -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _bpmcomplete BPM_complete.f90
```
(<\your\path\to\f2py.py>: most likely C:\Python27\Scripts\f2py.py)

## Running the Python code

This python code can be run from another file using:
```python
from BPM_complete_bladeseg import turbinepos
turbinepos(x,y,obs,wind,rpm,L,windvel,B,h,Rhub,rad,c,alpha,corr,py_for) # SPL at specified observer location (xo,yo,zo) for a given turbine rotation rate, blade geometry, and free stream wind speed
```
