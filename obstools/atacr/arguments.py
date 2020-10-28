# Copyright 2019 Pascal Audet & Helen Janiszewski
#
# This file is part of OBStools.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

Module containing the main utility functions used in the `OBStools` scripts
that accompany this package.

"""

# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_daylong_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_download_data.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used " +
        "to download and pre-process up to four-component " +
        "(H1, H2, Z and P), day-long seismograms to use in " +
        "noise corrections of vertical component of OBS data. " +
        "Data are requested from the internet using the client " +
        "services framework for a given date range. The stations " +
        "are processed one by one and the data are stored to disk.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma-separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the dictionary. " +
        "For instance, providing IU will match with all stations " +
        "in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_argument(
        "-C", "--channels",
        action="store",
        type=str,
        dest="channels",
        default="",
        help="Specify a comma-separated list of channels for " +
        "which to perform the transfer function analysis. " +
        "Possible options are H (for horizontal channels) or P " +
        "(for pressure channel). Specifying H allows " +
        "for tilt correction. Specifying P allows for compliance " +
        "correction. [Default looks for both horizontal and " +
        "pressure and allows for both tilt AND compliance corrections]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, " +
        "NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. " +
        "[Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to access and download " +
        "restricted data. [Default no user and password]")

    """
    # # Database Settings
    # DataGroup = parser.add_argument_group(parser, title="Local Data Settings", description="Settings associated with defining " \
    #     "and using a local data base of pre-downloaded day-long SAC files.")
    # DataGroup.add_argument("--local-data", action="store", type=str, dest="localdata", default=None, \
    #     help="Specify a comma separated list of paths containing day-long sac files of data already downloaded. " \
    #     "If data exists for a seismogram is already present on disk, it is selected preferentially over downloading " \
    #     "the data using the Client interface")
    # DataGroup.add_argument("--no-data-zero", action="store_true", dest="ndval", default=False, \
    #     help="Specify to force missing data to be set as zero, rather than default behaviour which sets to nan.")
    """

    # Constants Settings
    FreqGroup = parser.add_argument_group(
        title='Frequency Settings',
        description="Miscellaneous frequency settings")
    FreqGroup.add_argument(
        "--sampling-rate",
        action="store",
        type=float,
        dest="new_sampling_rate",
        default=5.,
        help="Specify new sampling rate (float, in Hz). [Default 5.]")
    FreqGroup.add_argument(
        "--units",
        action="store",
        type=str,
        dest="units",
        default="DISP",
        help="Choose the output seismogram units. Options are: " +
        "'DISP', 'VEL', 'ACC'. [Default 'DISP']")
    FreqGroup.add_argument(
        "--pre-filt",
        action="store",
        type=str,
        dest="pre_filt",
        default=None,
        help="Specify four comma-separated corner frequencies " +
        "(float, in Hz) for deconvolution pre-filter. " +
        "[Default 0.001,0.005,45.,50.]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with searching " +
        "for day-long seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start day for the data search. This will override any " +
        "station start times. " +
        "[Default start date for each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the event search. This will override any " +
        "station end times [Default end date for each station in database]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # create channel list
    if len(args.channels) > 0:
        args.channels = args.channels.split(',')
    else:
        args.channels = ["H", "P"]
    for cha in args.channels:
        if cha not in ["H", "P"]:
            parser.error("Error: Channel not recognized ", cha)

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for " +
                "User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # # Parse Local Data directories
    # if args.localdata is not None:
    #     args.localdata = args.localdata.split(',')
    # else:
    #     args.localdata = []

    # # Check NoData Value
    # if args.ndval:
    #     args.ndval = 0.0
    # else:
    #     args.ndval = nan

    if not type(args.new_sampling_rate) is float:
        raise(Exception(
            "Error: Type of --sampling-rate is not a float"))

    if args.units not in ['DISP', 'VEL', 'ACC']:
        raise(Exception(
            "Error: invalid --units argument. Choose among " +
            "'DISP', 'VEL', or 'ACC'"))
    if args.pre_filt is None:
        args.pre_filt = [0.001, 0.005, 45., 50.]
    else:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 4:
            raise(Exception(
                "Error: --pre-filt should contain 4 comma-separated floats"))

    return args


def get_event_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_download_event.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used " +
        "to download and pre-process four-component " +
        "(H1, H2, Z and P), two-hour-long seismograms for " +
        "individual events on which to apply the de-noising " +
        "algorithms. Data are requested from the internet using " +
        "the client services framework for a given date range. " +
        "The stations are processed one by one and the data are " +
        "stored to disk.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the "
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network [Default processes " +
        "all stations in the database]")
    parser.add_argument(
        "-C", "--channels",
        action="store",
        type=str,
        dest="channels",
        default="",
        help="Specify a comma-separated list of channels for " +
        "which to perform the transfer function analysis. " +
        "Possible options are H (for horizontal channels) or P " +
        "(for pressure channel). Specifying H allows " +
        "for tilt correction. Specifying P allows for compliance " +
        "correction. [Default looks for both horizontal and " +
        "pressure and allows for both tilt AND compliance corrections]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: BGR, " +
        "ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, NEIP, " +
        "NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. [Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to access and download " +
        "restricted data. [Default no user and password]")

    """
    # # Database Settings
    # DataGroup = parser.add_argument_group(parser, title="Local Data Settings", description="Settings associated with defining " \
    #     "and using a local data base of pre-downloaded day-long SAC files.")
    # DataGroup.add_argument("--local-data", action="store", type=str, dest="localdata", default=None, \
    #     help="Specify a comma separated list of paths containing day-long sac files of data already downloaded. " \
    #     "If data exists for a seismogram is already present on disk, it is selected preferentially over downloading " \
    #     "the data using the Client interface")
    # DataGroup.add_argument("--no-data-zero", action="store_true", dest="ndval", default=False, \
    #     help="Specify to force missing data to be set as zero, rather than default behaviour which sets to nan.")
    """

    # Constants Settings
    FreqGroup = parser.add_argument_group(
        title='Frequency Settings',
        description="Miscellaneous frequency settings")
    FreqGroup.add_argument(
        "--sampling-rate",
        action="store",
        type=float,
        dest="new_sampling_rate",
        default=5.,
        help="Specify new sampling rate (float, in Hz). " +
        "[Default 5.]")
    FreqGroup.add_argument(
        "--pre-filt",
        action="store",
        type=str,
        dest="pre_filt",
        default=None,
        help="Specify four comma-separated corner " +
        "frequencies (float, in Hz) for deconvolution " +
        "pre-filter. [Default 0.001,0.005,45.,50.]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with refining " +
        "the events to include in matching station " +
        "pairs")
    EventGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event " +
        "search. This will override any station start " +
        "times. [Default start date of each station in " +
        "database]")
    EventGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event " +
        "search. This will override any station end times " +
        "[Default end date of each station in database]")
    EventGroup.add_argument(
        "--reverse-order", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour " +
        "starts at oldest event and works towards most " +
        "recent. Specify reverse order and instead the " +
        "program will start with the most recent events " +
        "and work towards older")
    EventGroup.add_argument(
        "--min-mag",
        action="store",
        type=float,
        dest="minmag",
        default=5.5,
        help="Specify the minimum magnitude of event " +
        "for which to search. [Default 5.5]")
    EventGroup.add_argument(
        "--max-mag",
        action="store",
        type=float,
        dest="maxmag",
        default=None,
        help="Specify the maximum magnitude of event " +
        "for which to search. " +
        "[Default None, i.e. no limit]")

    # Geometry Settings
    GeomGroup = parser.add_argument_group(
        title="Geometry Settings",
        description="Settings associatd with the " +
        "event-station geometries")
    GeomGroup.add_argument(
        "--min-dist",
        action="store",
        type=float,
        dest="mindist",
        default=30.,
        help="Specify the minimum great circle distance " +
        "(degrees) between the station and event. " +
        "[Default 30]")
    GeomGroup.add_argument(
        "--max-dist",
        action="store",
        type=float,
        dest="maxdist",
        default=120.,
        help="Specify the maximum great circle distance " +
        "(degrees) between the station and event. " +
        "[Default 120]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # create channel list
    if len(args.channels) == 1:
        args.channels = args.stkeys.split(',')
    else:
        args.channels = ["H", "P"]
    for cha in args.channels:
        if cha not in ["H", "P"]:
            parser.error("Error: Channel not recognized ", cha)

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for User " +
                "Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # # Parse Local Data directories
    # if args.localdata is not None:
    #     args.localdata = args.localdata.split(',')
    # else:
    #     args.localdata = []

    # # Check NoData Value
    # if args.ndval:
    #     args.ndval = 0.0
    # else:
    #     args.ndval = nan

    if not type(args.new_sampling_rate) is float:
        raise(Exception("Error: Type of --sampling-rate is not a float"))

    if args.pre_filt is None:
        args.pre_filt = [0.001, 0.005, 45., 50.]
    else:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 4:
            raise(Exception(
                "Error: --pre-filt should contain 4 " +
                "comma-separated floats"))

    return args


def get_dailyspec_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_daily_spectra.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used "
        "to extract shorter windows from the day-long " +
        "seismograms, calculate the power-spectral " +
        "properties, flag windows for outlier PSDs and " +
        "calculate daily averages of the corresponding " +
        "Fourier transforms. The stations are processed " +
        "one by one and the data are stored to disk. The " +
        "program will look for data saved in the previous " +
        "steps and use all available components.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for day-long seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. " +
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. " +
        "This will override any station end times. " +
        "[Default end date of each station n database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_argument(
        "--window",
        action="store",
        type=float,
        dest="window",
        default=7200.,
        help="Specify window length in seconds. " +
        "Default value is highly recommended. "
        "Program may not be stable for large deviations " +
        "from default value. [Default 7200. (or 2 hours)]")
    ConstGroup.add_argument(
        "--overlap",
        action="store",
        type=float,
        dest="overlap",
        default=0.3,
        help="Specify fraction of overlap between windows. " +
        "[Default 0.3 (or 30%)]")
    ConstGroup.add_argument(
        "--minwin",
        action="store",
        type=int,
        dest="minwin",
        default=10,
        help="Specify minimum number of 'good' windows " +
        "in any given day to continue with analysis. " +
        "[Default 10]")
    ConstGroup.add_argument(
        "--freq-band",
        action="store",
        type=str,
        dest="pd",
        default=None,
        help="Specify comma-separated frequency limits " +
        "(float, in Hz) over which to calculate spectral " +
        "features used in flagging the bad windows. " +
        "[Default 0.004,2.0]")
    ConstGroup.add_argument(
        "--tolerance",
        action="store",
        type=float,
        dest="tol",
        default=2.0,
        help="Specify parameter for tolerance threshold. " +
        "If spectrum > std*tol, window is flagged as bad. " +
        "[Default 2.0]")
    ConstGroup.add_argument(
        "--alpha",
        action="store",
        type=float,
        dest="alpha",
        default=0.05,
        help="Specify confidence level for f-test, " +
        "for iterative flagging of windows. " +
        "[Default 0.05, or 95% confidence]")
    ConstGroup.add_argument(
        "--raw",
        action="store_true",
        dest="raw",
        default=False,
        help="Raw spectra will be used in calculating " +
        "spectral features for flagging. " +
        "[Default uses smoothed spectra]")
    ConstGroup.add_argument(
        "--no-rotation",
        action="store_false",
        dest="calc_rotation",
        default=True,
        help="Do not rotate horizontal components " +
        "to tilt direction. [Default calculates rotation]")

    # Constants Settings
    FigureGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_argument(
        "--figQC",
        action="store_true",
        dest="fig_QC",
        default=False,
        help="Plot Quality-Control figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--debug",
        action="store_true",
        dest="debug",
        default=False,
        help="Plot intermediate steps for debugging. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figAverage",
        action="store_true",
        dest="fig_average",
        default=False,
        help="Plot daily average figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figCoh",
        action="store_true",
        dest="fig_coh_ph",
        default=False,
        help="Plot Coherence and Phase figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--save-fig",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure(s). [Default " +
        "does not save figure]")
    FigureGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    if args.raw:
        args.smooth = False
    else:
        args.smooth = True

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Check input frequency band
    if args.pd is None:
        args.pd = [0.004, 2.0]
    else:
        args.pd = [float(val) for val in args.pd.split(',')]
        args.pd = sorted(args.pd)
        if (len(args.pd)) != 2:
            raise(Exception(
                "Error: --freq-band should contain 2 " +
                "comma-separated floats"))

    return args


def get_cleanspec_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_clean_spectra.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used "
        "to extract daily spectra calculated from " +
        "`obs_daily_spectra.py` and flag days for outlier " +
        "PSDs and calculate spectral averages of the " +
        "corresponding Fourier transforms over the entire " +
        "time period specified. The stations are processed " +
        "one by one and the data are stored to disk.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must " +
        "be contained within the station database. Partial " +
        "keys will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for day-long seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. " +
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. " +
        "This will override any station end times. " +
        "[Default end date of each station in database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_argument(
        "--freq-band",
        action="store",
        type=str,
        dest="pd",
        default=None,
        help="Specify comma-separated frequency limits " +
        "(float, in Hz) over which to calculate spectral " +
        "features used in flagging the days/windows. " +
        "[Default 0.004,2.0]")
    ConstGroup.add_argument(
        "--tolerance",
        action="store",
        type=float,
        dest="tol",
        default=1.5,
        help="Specify parameter for tolerance threshold. " +
        "If spectrum > std*tol, window is flagged as bad. " +
        "[Default 1.5]")
    ConstGroup.add_argument(
        "--alpha",
        action="store",
        type=float,
        dest="alpha",
        default=0.05,
        help="Confidence level for f-test, for iterative " +
        r"flagging of windows. [Default 0.05, or 95% confidence]")

    # Constants Settings
    FigureGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_argument(
        "--figQC",
        action="store_true",
        dest="fig_QC",
        default=False,
        help="Plot Quality-Control figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--debug",
        action="store_true",
        dest="debug",
        default=False,
        help="Plot intermediate steps for debugging. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figAverage",
        action="store_true",
        dest="fig_average",
        default=False,
        help="Plot daily average figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figCoh",
        action="store_true",
        dest="fig_coh_ph",
        default=False,
        help="Plot Coherence and Phase figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figCross",
        action="store_true",
        dest="fig_av_cross",
        default=False,
        help="Plot cross-spectra figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--save-fig",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure(s). [Default " +
        "does not save figure]")
    FigureGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    if args.pd is None:
        args.pd = [0.004, 2.0]
    else:
        args.pd = [float(val) for val in args.pd.split(',')]
        args.pd = sorted(args.pd)
        if (len(args.pd)) != 2:
            raise(Exception(
                "Error: --freq-band should contain 2 " +
                "comma-separated floats"))

    return args


def get_transfer_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_transfer functions.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used "
        "to calculate transfer functions between various " +
        "components, to be used in cleaning vertical " +
        "component of OBS data. The noise data can be " +
        "those obtained from the daily spectra (i.e., " +
        "from `obs_daily_spectra.py`) or those obtained " +
        "from the averaged noise spectra (i.e., from " +
        "`obs_clean_spectra.py`). Flags are available " +
        "to specify the source of data to use as well as " +
        "the time range over which to calculate the " +
        "transfer functions. The stations are processed " +
        "one by one and the data are stored to disk.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must " +
        "be contained within the station database. Partial " +
        "keys will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with searching "
        "for day-long seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. "
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. "
        "This will override any station end times. " +
        "[Default end date of each station in database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_argument(
        "--skip-daily",
        action="store_true",
        dest="skip_daily",
        default=False,
        help="Skip daily spectral averages in construction " +
        "of transfer functions. [Default False]")
    ConstGroup.add_argument(
        "--skip-clean",
        action="store_true",
        dest="skip_clean",
        default=False,
        help="Skip cleaned spectral averages in " +
        "construction of transfer functions. Defaults " +
        "to True if data cannot be found in default " +
        "directory. [Default False]")

    # Constants Settings
    FigureGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_argument(
        "--figTF",
        action="store_true",
        dest="fig_TF",
        default=False,
        help="Plot transfer function figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--save-fig",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure(s). [Default " +
        "does not save figure]")
    FigureGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "start time: " + args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "end time: " + args.endT)
    else:
        args.endT = None

    if args.skip_clean and args.skip_daily:
        parser.error(
            "Error: cannot skip both daily and clean averages")

    return args


def get_correct_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_correct_event.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used "
        "to extract transfer functions between various " +
        "components, and use them to clean vertical " +
        "component of OBS data for selected events. The " +
        "noise data can be those obtained from the daily " +
        "spectra (i.e., from `obs_daily_spectra.py`) "
        "or those obtained from the averaged noise spectra " +
        "(i.e., from `obs_clean_spectra.py`). Flags are " +
        "available to specify the source of data to use as " +
        "well as the time range for given events. "
        "The stations are processed one by one and the " +
        "data are stored to disk.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must be "
        "contained within the station database. Partial keys " +
        "will be used to match against those in the "
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network. [Default processes "
        "all stations in the database]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for specific event-related seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the event search. "
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event search. "
        "This will override any station end times. [Default " +
        "end date of each station in database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default " +
        "values and settings")
    ConstGroup.add_argument(
        "--skip-daily",
        action="store_true",
        dest="skip_daily",
        default=False,
        help="Skip daily spectral averages in application " +
        "of transfer functions. [Default False]")
    ConstGroup.add_argument(
        "--skip-clean",
        action="store_true",
        dest="skip_clean",
        default=False,
        help="Skip cleaned spectral averages in " +
        "application of transfer functions. " +
        "[Default False]")
    ConstGroup.add_argument(
        "--fmin",
        action="store",
        type=float,
        dest="fmin",
        default="0.006666666666666667",
        help="Low frequency corner (in Hz) for " +
        "plotting the raw (un-corrected) seismograms. "
        "Filter is a 2nd order, zero phase butterworth " +
        "filter. [Default 1./150.]")
    ConstGroup.add_argument(
        "--fmax",
        action="store",
        type=float,
        dest="fmax",
        default="0.1",
        help="High frequency corner (in Hz) for " +
        "plotting the raw (un-corrected) seismograms. "
        "Filter is a 2nd order, zero phase butterworth " +
        "filter. [Default 1./10.]")

    # Constants Settings
    FigureGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_argument(
        "--figRaw",
        action="store_true",
        dest="fig_event_raw",
        default=False,
        help="Plot raw seismogram figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figClean",
        action="store_true",
        dest="fig_plot_corrected",
        default=False,
        help="Plot cleaned vertical seismogram figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--save-fig",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure(s). [Default " +
        "does not save figure]")
    FigureGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "start time: " + args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "end time: " + args.endT)
    else:
        args.endT = None

    if args.skip_clean and args.skip_daily:
        parser.error(
            "Error: cannot skip both daily and clean averages")

    return args

