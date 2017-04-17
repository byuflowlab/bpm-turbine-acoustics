#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup, find_packages

setup(
    name='bpmacoustic',
    version='1.2.0',
    description='BPM acoustic model for turbine noise prediction',
    author='Eric B. Tingey',
    author_email='ebtingey@byu.edu',
    url='https://github.com/byuflowlab/bpm-turbine-acoustics',
    package_dir={'': 'acoustic_model'},
    py_modules=['bpmacoustic'],
    license='MIT License',
    zip_safe=False
)

from numpy.distutils.core import setup, Extension
setup(
    name='bpmacoustic',
    version='1.2.0',
    package_dir={'': 'acoustic_model'},
    ext_modules=[Extension('_bpmacoustic', ['acoustic_model/BPM_Acoustic_Model.f90'], extra_compile_args=['-O2'])],
)

setup(
    name='bpmvawtacoustic',
    version='1.3.0',
    package_dir={'': 'acoustic_model'},
    ext_modules=[Extension('_bpmvawtacoustic', ['acoustic_model/BPM_VAWT_Acoustic_Model.f90'], extra_compile_args=['-O2'])],
)
