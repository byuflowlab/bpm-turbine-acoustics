# BPM Turbine Acoustics

Turbine acoustic code using the BPM equations developed by Brooks, Pope, and Marcolini

Developed by Eric Tingey at Brigham Young University, 2015

This code models the acoustic propagation of a wind turbine based on turbulent boundary layer edge noise, separation stall noise, tip vortex formation noise, laminar boundary layer vortex shedding noise, and trailing edge bluntness vortex shedding noise. Turbulent inflow noise is not assumed in this current code. The semi-empiracal equations were developed from the NACA 0012 airfoil data and the blade segments used in the test file are based on the NREL 5-MW wind turbine. Scaling of the segments is based on the blade length specified.

Brooks, T., Pope, D., and Marcolini, M., “Aipower Self-Noise and Prediction,” NASA, 1989.


## Installation instructions

- system requirements: gfortran (using MinGW for Windows in order to use the commands here), python 2.7, numpy, scipy
- run:
```
python setup.py install
```
- or navigate to the directory and run the following command in the terminal to build the Fortran code:

Mac
```
$ cd wake_model
$ f2py -c  --opt=-O2 -m _bpmacoustic BPM_Acoustic_Model.f90
```

Windows
```
cd wake_model
python <\your\path\to\f2py.py> -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _bpmacoustic BPM_Acoustic_Model.f90
```
(<\your\path\to\f2py.py>: most likely C:\Python27\Scripts\f2py.py)

## Running the Python code

This python code can be run from another file using:
```python
SPL = _bpmacoustic.turbinepos(x, y, obs, winddir, rpm, windvel, B, h, rad, c, alpha, noise_corr) # SPL at specified observer location (xo,yo,zo) for a given turbine rotation rate, blade geometry, and free stream wind speed
```
