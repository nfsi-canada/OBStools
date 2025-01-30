import setuptools
import os.path
from os import listdir
import re
from setuptools import setup
from pathlib import Path


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname, encoding='utf-8') as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


# scripts = [str(x) for x in Path('Scripts').iterdir() if x.is_file()]

setup(
    name='obstools',
    version=find_version('obstools', '__init__.py'),
    description='Python tools for ocean bottom seismic instruments',
    author='Pascal Audet, Helen Janiszewski',
    author_email='pascal.audet@uottawa.ca',
    maintainer='Pascal Audet, Helen Janiszewski',
    maintainer_email='pascal.audet@uottawa.ca, hajanisz@hawaii.edu',
    url='https://github.com/nfsi-canada/OBStools/archive/OBStools-0.1.4.tar.gz',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'],
    install_requires=['numpy', 'obspy', 'stdb', 'scipy', 'pandas'],
    python_requires='>=3.6',
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'': ['examples/meta/*.pkl', 'examples/event/*.pkl',
                       'examples/data/*.SAC', 'examples/event/*.SAC']},
    #    scripts=scripts)
    entry_points={
        'console_scripts':
        ['atacr_download_data=obstools.scripts.atacr_download_data:main',
         # 'atacr_download_data_xml=obstools.scripts.atacr_download_data_xml:main',
         'atacr_download_event=obstools.scripts.atacr_download_event:main',
         'atacr_daily_spectra=obstools.scripts.atacr_daily_spectra:main',
         'atacr_clean_spectra=obstools.scripts.atacr_clean_spectra:main',
         'atacr_transfer_functions=obstools.scripts.atacr_transfer_functions:main',
         'atacr_correct_event=obstools.scripts.atacr_correct_event:main',
         'comply_calculate=obstools.scripts.comply_calculate:main']})
