# BPM Acoustics Optimization Code

This folder contains the optimization codes used to run the BPM acoustics model with the FLORIS turbine wake model. After building the codes (see below), run the optimization from the command line using:

    $ python Optimization.py [SPL limit value (dB)]

The optimization can also be run using a hardcoded SPL limit found at line 271 and commenting out line 270. Two example wind farms (Lissett Airfield Wind Farm and Rosiere Wind Farm) can be used as the initial layouts can can be chosen at lines 274 and 275. For faster SPL calculations, the joblib module is implemented and can be commented in or out on lines 120 and 121 (https://pythonhosted.org/joblib/).

This older version of the FLORIS model was used in the optimizations used for our specific publication, but a new version of the codes along with more detailed instructions at: https://github.com/wisdem/florisse. Also required for the optimization is CCBlade, of which the older version using in our optimizations, is included, but the newer version can be found at: https://github.com/wisdem/ccblade. All optimizations were run with SNOPT, a Python optimization code, which is also required to run the optimization code. The following instructions outline how to build the FLORIS codes for use in reference with the BPM code (copied from https://github.com/WISDEM/FLORISSE/tree/develop). The newer adjustments require joblib to run the code with parallelization.

Implementation of the FLOw Redirection and Induction in Steady-state (FLORIS) wind farm wake model for WISDEM, written in Python.

Created by Pieter Gebraad and Paul Fleming. Copyright (c) NREL. All rights reserved.

REQUIRED PYTHON LIBRARIES:
- OpenMDAO
- Numpy
- WISDEM/akima

-For a description of the FLORIS model refer to: P. M. O. Gebraad, F. W. Teeuwisse, J.-W. van Wingerden, P. A. Fleming, S. D. Ruben, J. R. Marden, and L. Pao, “Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study,” in Wind Energy, 2014, DOI: 10.1002/we.1822.

## Installation instructions MAC  
### system requirements  
    gfortran  
    gcc  
    python 2.7.x  
    numpy  
    openmdao >= v1.5  
### run the following commands from src/florisse:  
    $ gfortran -c adBuffer.f  
    $ gcc -c adStack.c  
    $ f2py -c --opt=-O2 -m _florisunify florisUnified.f90 adBuffer.o adStack.o  


## Installation instructions Windows  
### system requirements  
    gfortran  
    gcc  
    mingw  
    python 2.7.x  
    numpy  
    openmdao >= v01.5
### run the following commands from src\florisse:  
    $ gfortran -c adBuffer.f  
    $ gcc -c adStack.c  
    $ python \your\path\to\f2py.py -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _florisunify florisUnified.f90 adBuffer.o adStack.o  
        (most likely your path is C:\python27\Scripts\f2py.py)  
### if you get an error in the line "as=b['args']" try to update numpy
    ($ pip install numpy --upgrade)  


## Installation instructions Marylou (BYU Supercomputer)
### module dependencies ($ module load <module name>)(purge all modules first)  
    petsc/3.6.3  
    python/2/7  
    compiler_gnu/4.9.2  
    mpi/openmpi-1.8.4_gnu-4.9.2  
### python dependencies ($ pip install --user <package name>)  
    openmdao >= v1.5 (use a clone of the repo and {$ pip install --user -e .} from top level of
              acquired repo)  
    mpi4py  
    petsc4py      
### compiler FLORIS (clone with ssh on Marylou)  
    $ cd src/florisse  
    $ gcc -fPIC -c adStack.c  
    $ gfortran -fPIC -c adBuffer.f  
    $ f2py -c --opt=-O2 -m _florisunify florisUnified.f90 adBuffer.o adStack.o   
