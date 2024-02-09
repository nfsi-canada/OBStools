#---I don't have a spot for these in the package yet. Still workging on the PSD testing:
# import scipy.stats as stats
# import scipy as sc
# from numpy.lib.stride_tricks import sliding_window_view

#---Not used in my QA. Just have it here for reference:
# from obstools.atacr.plotting import fig_QC, fig_average, fig_av_cross, fig_coh_ph, fig_TF, fig_comply, fig_event_raw, fig_event_corrected

# -------Import main-------::
import ObsQA

# -------Imports for plot code-------::
import matplotlib.pyplot as plt
import obspy
from obspy import taup
import re
from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
import numpy as np
import folium
from branca.element import Figure

# -------Imports for io code-------::
import obstools
from obstools.atacr import DayNoise, TFNoise, EventStream, StaNoise, utils
import obstools.atacr.plotting as atplot
from obstools.scripts import comply_calculate, atacr_clean_spectra, atacr_correct_event, atacr_daily_spectra, atacr_download_data, atacr_download_event, atacr_transfer_functions
from stdb.scripts import query_fdsn_stdb
from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
import scipy.io as spio
from scipy.signal import csd
import glob as g
import pandas as pd
import numpy as np
import pickle as pkl
import pickle
from obspy.clients.fdsn import Client
import datetime
import os
import sys

import obspy

import geemap
import ee
import matplotlib
import distinctipy
from branca.element import Template, MacroElement

import scipy.io as spio
import numpy as np
import pandas as pd
import glob as g