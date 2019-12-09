import os.path
from os import listdir
import re
from numpy.distutils.core import setup

def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")

scripts = ['Scripts/' + i for i in listdir('Scripts/')]

setup(
	name='obstools',
    version=find_version('obstools', '__init__.py'),
	description='Python tools for ocean bottom seismic instruments',
	author='Pascal Audet, Helen Janiszewski',
    author_email='pascal.audet@uottawa.ca',
	maintainer='Pascal Audet, Helen Janiszewski',
	maintainer_email='pascal.audet@uottawa.ca, hajanisz@hawaii.edu',
    url='https://github.com/paudetseis/OBStools',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    install_requires=['numpy', 'obspy', 'stdb'],
    python_requires='>=3.6',
	packages=['obstools'],
	scripts=scripts)
