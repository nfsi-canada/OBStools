#!/usr/bin/env python

from distutils.core import setup
from os import listdir

scripts = ['Scripts/' + i for i in listdir('Scripts/')]

setup(
	name =             'obstools',
	version =          '0.0.1',
	description =      'Python Module for managing Station Databases',
	author =           'Andrew Schaeffer, Pascal Audet',
	maintainer =       'Andrew Schaeffer, Pascal Audet',
	maintainer_email = 'andrew.schaeffer@canada.ca, pascal.audet@uottawa.ca',
    classifiers      = [
                         'Development Status :: 3 - Alpha',
                         'License :: OSI Approved :: MIT License',
                         'Programming Language :: Python :: 3.6',
                         'Programming Language :: Python :: 3.7'
                       ],
    install_requires = ['obspy', 'stdb'],
    python_requires =  '>=3.6',
	packages =         ['obstools'],
	scripts =          scripts,
	url =              '')
