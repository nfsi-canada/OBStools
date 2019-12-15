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
from obspy import UTCDateTime
from numpy import nan, isnan


def get_daylong_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_download_data.py` that accompany this package.

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used " +
        "to download and pre-process up to four-component " +
        "(H1, H2, Z and P), day-long seismograms to use in " +
        "noise corrections of vertical component of OBS data. " +
        "Data are requested from the internet using the client " +
        "services framework for a given date range. The stations " +
        "are processed one by one and the data are stored to disk.")

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type="string",
        dest="stkeys",
        default="",
        help="Specify a comma-separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the dictionary. " +
        "For instance, providing IU will match with all stations " +
        "in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_option(
        "-C", "--channels",
        action="store",
        type="string",
        dest="channels",
        default="",
        help="Specify a comma-separated list of channels for " +
        "which to perform the transfer function analysis. " +
        "Possible options are H (for horizontal channels) or P " +
        "(for pressure channel). Specifying H allows " +
        "for tilt correction. Specifying P allows for compliance " +
        "correction. [Default looks for both horizontal and " +
        "pressure and allows for both tilt AND compliance corrections]")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Server Settings
    ServerGroup = OptionGroup(
        parser,
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_option(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, " +
        "NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. " +
        "[Default IRIS]")
    ServerGroup.add_option(
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
    # DataGroup = OptionGroup(parser, title="Local Data Settings", description="Settings associated with defining " \
    #     "and using a local data base of pre-downloaded day-long SAC files.")
    # DataGroup.add_option("--local-data", action="store", type="string", dest="localdata", default=None, \
    #     help="Specify a comma separated list of paths containing day-long sac files of data already downloaded. " \
    #     "If data exists for a seismogram is already present on disk, it is selected preferentially over downloading " \
    #     "the data using the Client interface")
    # DataGroup.add_option("--no-data-zero", action="store_true", dest="ndval", default=False, \
    #     help="Specify to force missing data to be set as zero, rather than default behaviour which sets to nan.")
    """

    # Constants Settings
    FreqGroup = OptionGroup(
        parser, title='Frequency Settings',
        description="Miscellaneous frequency settings")
    FreqGroup.add_option(
        "--sampling-rate",
        action="store",
        type="float",
        dest="new_sampling_rate",
        default=5.,
        help="Specify new sampling rate (float, in Hz). [Default 5.]")
    FreqGroup.add_option(
        "--pre-filt",
        action="store",
        type="string",
        dest="pre_filt",
        default=None,
        help="Specify four comma-separated corner frequencies " +
        "(float, in Hz) for deconvolution pre-filter. " +
        "[Default 0.001,0.005,45.,50.]")

    # Event Selection Criteria
    DaysGroup = OptionGroup(
        parser,
        title="Time Search Settings",
        description="Time settings associated with searching " +
        "for day-long seismograms")
    DaysGroup.add_option(
        "--start",
        action="store",
        type="string",
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start day for the data search. This will override any " +
        "station start times. " +
        "[Default start date for each station in database]")
    DaysGroup.add_option(
        "--end",
        action="store",
        type="string",
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the event search. This will override any " +
        "station end times [Default end date for each station in database]")

    parser.add_option_group(ServerGroup)
    # parser.add_option_group(DataGroup)
    parser.add_option_group(DaysGroup)
    parser.add_option_group(FreqGroup)
    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Error: Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Error: Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    # create channel list
    if len(opts.channels) > 0:
        opts.channels = opts.channels.split(',')
    else:
        opts.channels = ["H", "P"]
    for cha in opts.channels:
        if cha not in ["H", "P"]:
            parser.error("Error: Channel not recognized ", cha)

    # construct start time
    if len(opts.startT) > 0:
        try:
            opts.startT = UTCDateTime(opts.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                opts.startT)
    else:
        opts.startT = None

    # construct end time
    if len(opts.endT) > 0:
        try:
            opts.endT = UTCDateTime(opts.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                opts.endT)
    else:
        opts.endT = None

    # Parse User Authentification
    if not len(opts.UserAuth) == 0:
        tt = opts.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for " +
                "User Authentification")
        else:
            opts.UserAuth = tt
    else:
        opts.UserAuth = []

    # # Parse Local Data directories
    # if opts.localdata is not None:
    #     opts.localdata = opts.localdata.split(',')
    # else:
    #     opts.localdata = []

    # # Check NoData Value
    # if opts.ndval:
    #     opts.ndval = 0.0
    # else:
    #     opts.ndval = nan

    if not type(opts.new_sampling_rate) is float:
        raise(Exception(
            "Error: Type of --sampling-rate is not a float"))

    if opts.pre_filt is None:
        opts.pre_filt = [0.001, 0.005, 45., 50.]
    else:
        opts.pre_filt = [float(opts.pre_filt.split(','))]
        opts.pre_filt = sorted(opts.pre_filt)
        if (len(opts.pre_filt)) != 4:
            raise(Exception(
                "Error: --pre-filt should contain 4 comma-separated floats"))

    return (opts, indb)


def get_event_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_download_event.py` that accompany this package.

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used " +
        "to download and pre-process four-component " +
        "(H1, H2, Z and P), two-hour-long seismograms for " +
        "individual events on which to apply the de-noising " +
        "algorithms. Data are requested from the internet using " +
        "the client services framework for a given date range. " +
        "The stations are processed one by one and the data are " +
        "stored to disk.")

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type="string",
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the "
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network [Default processes " +
        "all stations in the database]")
    parser.add_option(
        "-C", "--channels",
        action="store",
        type="string",
        dest="channels",
        default="",
        help="Specify a comma-separated list of channels for " +
        "which to perform the transfer function analysis. " +
        "Possible options are H (for horizontal channels) or P " +
        "(for pressure channel). Specifying H allows " +
        "for tilt correction. Specifying P allows for compliance " +
        "correction. [Default looks for both horizontal and " +
        "pressure and allows for both tilt AND compliance corrections]")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Server Settings
    ServerGroup = OptionGroup(
        parser,
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_option(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: BGR, " +
        "ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, NEIP, " +
        "NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. [Default IRIS]")
    ServerGroup.add_option(
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
    # DataGroup = OptionGroup(parser, title="Local Data Settings", description="Settings associated with defining " \
    #     "and using a local data base of pre-downloaded day-long SAC files.")
    # DataGroup.add_option("--local-data", action="store", type="string", dest="localdata", default=None, \
    #     help="Specify a comma separated list of paths containing day-long sac files of data already downloaded. " \
    #     "If data exists for a seismogram is already present on disk, it is selected preferentially over downloading " \
    #     "the data using the Client interface")
    # DataGroup.add_option("--no-data-zero", action="store_true", dest="ndval", default=False, \
    #     help="Specify to force missing data to be set as zero, rather than default behaviour which sets to nan.")
    """

    # Constants Settings
    FreqGroup = OptionGroup(
        parser,
        title='Frequency Settings',
        description="Miscellaneous frequency settings")
    FreqGroup.add_option(
        "--sampling-rate",
        action="store",
        type="float",
        dest="new_sampling_rate",
        default=5.,
        help="Specify new sampling rate (float, in Hz). " +
        "[Default 5.]")
    FreqGroup.add_option(
        "--pre-filt",
        action="store",
        type="string",
        dest="pre_filt",
        default=None,
        help="Specify four comma-separated corner " +
        "frequencies (float, in Hz) for deconvolution " +
        "pre-filter. [Default 0.001,0.005,45.,50.]")

    # Event Selection Criteria
    EventGroup = OptionGroup(
        parser,
        title="Event Settings",
        description="Settings associated with refining " +
        "the events to include in matching station " +
        "pairs")
    EventGroup.add_option(
        "--start",
        action="store",
        type="string",
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event " +
        "search. This will override any station start " +
        "times. [Default start date of each station in " +
        "database]")
    EventGroup.add_option(
        "--end",
        action="store",
        type="string",
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event " +
        "search. This will override any station end times " +
        "[Default end date of each station in database]")
    EventGroup.add_option(
        "--reverse-order", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour " +
        "starts at oldest event and works towards most " +
        "recent. Specify reverse order and instead the " +
        "program will start with the most recent events " +
        "and work towards older")
    EventGroup.add_option(
        "--min-mag",
        action="store",
        type="float",
        dest="minmag",
        default=5.5,
        help="Specify the minimum magnitude of event " +
        "for which to search. [Default 5.5]")
    EventGroup.add_option(
        "--max-mag",
        action="store",
        type="float",
        dest="maxmag",
        default=None,
        help="Specify the maximum magnitude of event " +
        "for which to search. " +
        "[Default None, i.e. no limit]")

    # Geometry Settings
    GeomGroup = OptionGroup(
        parser,
        title="Geometry Settings",
        description="Settings associatd with the " +
        "event-station geometries")
    GeomGroup.add_option(
        "--min-dist",
        action="store",
        type="float",
        dest="mindist",
        default=30.,
        help="Specify the minimum great circle distance " +
        "(degrees) between the station and event. " +
        "[Default 30]")
    GeomGroup.add_option(
        "--max-dist",
        action="store",
        type="float",
        dest="maxdist",
        default=120.,
        help="Specify the maximum great circle distance " +
        "(degrees) between the station and event. " +
        "[Default 120]")

    parser.add_option_group(ServerGroup)
    # parser.add_option_group(DataGroup)
    parser.add_option_group(EventGroup)
    parser.add_option_group(GeomGroup)
    parser.add_option_group(FreqGroup)
    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Error: Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Error: Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    # create channel list
    if len(opts.channels) == 1:
        opts.channels = opts.stkeys.split(',')
    else:
        opts.channels = ["H", "P"]
    for cha in opts.channels:
        if cha not in ["H", "P"]:
            parser.error("Error: Channel not recognized ", cha)

    # construct start time
    if len(opts.startT) > 0:
        try:
            opts.startT = UTCDateTime(opts.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                opts.startT)
    else:
        opts.startT = None

    # construct end time
    if len(opts.endT) > 0:
        try:
            opts.endT = UTCDateTime(opts.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                opts.endT)
    else:
        opts.endT = None

    # Parse User Authentification
    if not len(opts.UserAuth) == 0:
        tt = opts.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for User " +
                "Authentification")
        else:
            opts.UserAuth = tt
    else:
        opts.UserAuth = []

    # # Parse Local Data directories
    # if opts.localdata is not None:
    #     opts.localdata = opts.localdata.split(',')
    # else:
    #     opts.localdata = []

    # # Check NoData Value
    # if opts.ndval:
    #     opts.ndval = 0.0
    # else:
    #     opts.ndval = nan

    if not type(opts.new_sampling_rate) is float:
        raise(Exception("Error: Type of --sampling-rate is not a float"))

    if opts.pre_filt is None:
        opts.pre_filt = [0.001, 0.005, 45., 50.]
    else:
        opts.pre_filt = [float(opts.pre_filt.split(','))]
        opts.pre_filt = sorted(opts.pre_filt)
        if (len(opts.pre_filt)) != 4:
            raise(Exception(
                "Error: --pre-filt should contain 4 " +
                "comma-separated floats"))

    return (opts, indb)


def get_dailyspec_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_daily_spectra.py` that accompany this package.

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used "
        "to extract shorter windows from the day-long " +
        "seismograms, calculate the power-spectral " +
        "properties, flag windows for outlier PSDs and " +
        "calculate daily averages of the corresponding " +
        "Fourier transforms. The stations are processed " +
        "one by one and the data are stored to disk. The " +
        "program will look for data saved in the previous " +
        "steps and use all available components.")

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type="string",
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = OptionGroup(
        parser,
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for day-long seismograms")
    DaysGroup.add_option(
        "--start",
        action="store",
        type="string",
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. " +
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_option(
        "--end",
        action="store",
        type="string",
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. " +
        "This will override any station end times. " +
        "[Default end date of each station n database]")

    # Constants Settings
    ConstGroup = OptionGroup(
        parser, title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_option(
        "--window",
        action="store",
        type="float",
        dest="window",
        default=7200.,
        help="Specify window length in seconds. " +
        "Default value is highly recommended. "
        "Program may not be stable for large deviations " +
        "from default value. [Default 7200. (or 2 hours)]")
    ConstGroup.add_option(
        "--overlap",
        action="store",
        type="float",
        dest="overlap",
        default=0.3,
        help="Specify fraction of overlap between windows. " +
        "[Default 0.3 (or 30%)]")
    ConstGroup.add_option(
        "--minwin",
        action="store",
        type="int",
        dest="minwin",
        default=10,
        help="Specify minimum number of 'good' windows " +
        "in any given day to continue with analysis. " +
        "[Default 10]")
    ConstGroup.add_option(
        "--freq-band",
        action="store",
        type="string",
        dest="pd",
        default=None,
        help="Specify comma-separated frequency limits " +
        "(float, in Hz) over which to calculate spectral " +
        "features used in flagging the bad windows. " +
        "[Default 0.004,2.0]")
    ConstGroup.add_option(
        "--tolerance",
        action="store",
        type="float",
        dest="tol",
        default=2.0,
        help="Specify parameter for tolerance threshold. " +
        "If spectrum > std*tol, window is flagged as bad. " +
        "[Default 2.0]")
    ConstGroup.add_option(
        "--alpha",
        action="store",
        type="float",
        dest="alpha",
        default=0.05,
        help="Specify confidence level for f-test, " +
        "for iterative flagging of windows. " +
        "[Default 0.05, or 95% confidence]")
    ConstGroup.add_option(
        "--raw",
        action="store_true",
        dest="raw",
        default=False,
        help="Raw spectra will be used in calculating " +
        "spectral features for flagging. " +
        "[Default uses smoothed spectra]")
    ConstGroup.add_option(
        "--no-rotation",
        action="store_false",
        dest="calc_rotation",
        default=True,
        help="Do not rotate horizontal components " +
        "to tilt direction. [Default calculates rotation]")

    # Constants Settings
    FigureGroup = OptionGroup(
        parser,
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_option(
        "--figQC",
        action="store_true",
        dest="fig_QC",
        default=False,
        help="Plot Quality-Control figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--debug",
        action="store_true",
        dest="debug",
        default=False,
        help="Plot intermediate steps for debugging. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--figAverage",
        action="store_true",
        dest="fig_average",
        default=False,
        help="Plot daily average figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--figCoh",
        action="store_true",
        dest="fig_coh_ph",
        default=False,
        help="Plot Coherence and Phase figure. " +
        "[Default does not plot figure]")

    parser.add_option_group(ConstGroup)
    parser.add_option_group(FigureGroup)
    parser.add_option_group(DaysGroup)
    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Error: Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Error: Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    if opts.raw:
        opts.smooth = False
    else:
        opts.smooth = True

    # construct start time
    if len(opts.startT) > 0:
        try:
            opts.startT = UTCDateTime(opts.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                opts.startT)
    else:
        opts.startT = None

    # construct end time
    if len(opts.endT) > 0:
        try:
            opts.endT = UTCDateTime(opts.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                opts.endT)
    else:
        opts.endT = None

    # Check input frequency band
    if opts.pd is None:
        opts.pd = [0.004, 2.0]
    else:
        opts.pd = [float(opts.pd.split(','))]
        opts.pd = sorted(opts.pd)
        if (len(opts.pd)) != 2:
            raise(Exception(
                "Error: --freq-band should contain 2 " +
                "comma-separated floats"))

    return (opts, indb)


def get_cleanspec_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_clean_spectra.py` that accompany this package.

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used "
        "to extract daily spectra calculated from " +
        "`obs_daily_spectra.py` and flag days for outlier " +
        "PSDs and calculate spectral averages of the " +
        "corresponding Fourier transforms over the entire " +
        "time period specified. The stations are processed " +
        "one by one and the data are stored to disk.")

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type="string",
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must " +
        "be contained within the station database. Partial " +
        "keys will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = OptionGroup(
        parser,
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for day-long seismograms")
    DaysGroup.add_option(
        "--start",
        action="store",
        type="string",
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. " +
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_option(
        "--end",
        action="store",
        type="string",
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. " +
        "This will override any station end times. " +
        "[Default end date of each station in database]")

    # Constants Settings
    ConstGroup = OptionGroup(
        parser,
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_option(
        "--freq-band",
        action="store",
        type="string",
        dest="pd",
        default=None,
        help="Specify comma-separated frequency limits " +
        "(float, in Hz) over which to calculate spectral " +
        "features used in flagging the days/windows. " +
        "[Default 0.004,2.0]")
    ConstGroup.add_option(
        "--tolerance",
        action="store",
        type="float",
        dest="tol",
        default=1.5,
        help="Specify parameter for tolerance threshold. " +
        "If spectrum > std*tol, window is flagged as bad. " +
        "[Default 1.5]")
    ConstGroup.add_option(
        "--alpha",
        action="store",
        type="float",
        dest="alpha",
        default=0.05,
        help="Confidence level for f-test, for iterative " +
        r"flagging of windows. [Default 0.05, or 95% confidence]")

    # Constants Settings
    FigureGroup = OptionGroup(
        parser, title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_option(
        "--figQC",
        action="store_true",
        dest="fig_QC",
        default=False,
        help="Plot Quality-Control figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--debug",
        action="store_true",
        dest="debug",
        default=False,
        help="Plot intermediate steps for debugging. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--figAverage",
        action="store_true",
        dest="fig_average",
        default=False,
        help="Plot daily average figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--figCoh",
        action="store_true",
        dest="fig_coh_ph",
        default=False,
        help="Plot Coherence and Phase figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--figCross",
        action="store_true",
        dest="fig_av_cross",
        default=False,
        help="Plot cross-spectra figure. " +
        "[Default does not plot figure]")

    parser.add_option_group(ConstGroup)
    parser.add_option_group(FigureGroup)
    parser.add_option_group(DaysGroup)
    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Error: Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Error: Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    # construct start time
    if len(opts.startT) > 0:
        try:
            opts.startT = UTCDateTime(opts.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                opts.startT)
    else:
        opts.startT = None

    # construct end time
    if len(opts.endT) > 0:
        try:
            opts.endT = UTCDateTime(opts.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                opts.endT)
    else:
        opts.endT = None

    if opts.pd is None:
        opts.pd = [0.004, 2.0]
    else:
        opts.pd = [float(opts.pd.split(','))]
        opts.pd = sorted(opts.pd)
        if (len(opts.pd)) != 2:
            raise(Exception(
                "Error: --freq-band should contain 2 " +
                "comma-separated floats"))

    return (opts, indb)


def get_transfer_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_transfer functions.py` that accompany this package.

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
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

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type="string",
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must " +
        "be contained within the station database. Partial " +
        "keys will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = OptionGroup(
        parser,
        title="Time Search Settings",
        description="Time settings associated with searching "
        "for day-long seismograms")
    DaysGroup.add_option(
        "--start",
        action="store",
        type="string",
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. "
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_option(
        "--end",
        action="store",
        type="string",
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. "
        "This will override any station end times. " +
        "[Default end date of each station in database]")

    # Constants Settings
    ConstGroup = OptionGroup(
        parser,
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_option(
        "--skip-daily",
        action="store_true",
        dest="skip_daily",
        default=False,
        help="Skip daily spectral averages in construction " +
        "of transfer functions. [Default False]")
    ConstGroup.add_option(
        "--skip-clean",
        action="store_true",
        dest="skip_clean",
        default=False,
        help="Skip cleaned spectral averages in " +
        "construction of transfer functions. Defaults " +
        "to True if data cannot be found in default " +
        "directory. [Default False]")

    # Constants Settings
    FigureGroup = OptionGroup(
        parser,
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_option(
        "--figTF",
        action="store_true",
        dest="fig_TF",
        default=False,
        help="Plot transfer function figure. " +
        "[Default does not plot figure]")

    parser.add_option_group(ConstGroup)
    parser.add_option_group(FigureGroup)
    parser.add_option_group(DaysGroup)
    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Error: Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Error: Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    # construct start time
    if len(opts.startT) > 0:
        try:
            opts.startT = UTCDateTime(opts.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "start time: " + opts.startT)
    else:
        opts.startT = None

    # construct end time
    if len(opts.endT) > 0:
        try:
            opts.endT = UTCDateTime(opts.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "end time: " + opts.endT)
    else:
        opts.endT = None

    if opts.skip_clean and opts.skip_daily:
        parser.error(
            "Error: cannot skip both daily and clean averages")

    return (opts, indb)


def get_correct_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_correct_event.py` that accompany this package.

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
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

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type="string",
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must be "
        "contained within the station database. Partial keys " +
        "will be used to match against those in the "
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network. [Default processes "
        "all stations in the database]")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = OptionGroup(
        parser,
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for specific event-related seismograms")
    DaysGroup.add_option(
        "--start",
        action="store",
        type="string",
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the event search. "
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_option(
        "--end",
        action="store",
        type="string",
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event search. "
        "This will override any station end times. [Default " +
        "end date of each station in database]")

    # Constants Settings
    ConstGroup = OptionGroup(
        parser,
        title='Parameter Settings',
        description="Miscellaneous default " +
        "values and settings")
    ConstGroup.add_option(
        "--skip-daily",
        action="store_true",
        dest="skip_daily",
        default=False,
        help="Skip daily spectral averages in application " +
        "of transfer functions. [Default False]")
    ConstGroup.add_option(
        "--skip-clean",
        action="store_true",
        dest="skip_clean",
        default=False,
        help="Skip cleaned spectral averages in " +
        "application of transfer functions. " +
        "[Default False]")
    ConstGroup.add_option(
        "--fmin",
        action="store",
        type="float",
        dest="fmin",
        default="0.006666666666666667",
        help="Low frequency corner (in Hz) for " +
        "plotting the raw (un-corrected) seismograms. "
        "Filter is a 2nd order, zero phase butterworth " +
        "filter. [Default 1./150.]")
    ConstGroup.add_option(
        "--fmax",
        action="store",
        type="float",
        dest="fmax",
        default="0.1",
        help="High frequency corner (in Hz) for " +
        "plotting the raw (un-corrected) seismograms. "
        "Filter is a 2nd order, zero phase butterworth " +
        "filter. [Default 1./10.]")

    # Constants Settings
    FigureGroup = OptionGroup(
        parser, title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_option(
        "--figRaw",
        action="store_true",
        dest="fig_event_raw",
        default=False,
        help="Plot raw seismogram figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_option(
        "--figClean",
        action="store_true",
        dest="fig_plot_corrected",
        default=False,
        help="Plot cleaned vertical seismogram figure. " +
        "[Default does not plot figure]")

    parser.add_option_group(ConstGroup)
    parser.add_option_group(FigureGroup)
    parser.add_option_group(DaysGroup)
    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Error: Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Error: Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    # construct start time
    if len(opts.startT) > 0:
        try:
            opts.startT = UTCDateTime(opts.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "start time: " + opts.startT)
    else:
        opts.startT = None

    # construct end time
    if len(opts.endT) > 0:
        try:
            opts.endT = UTCDateTime(opts.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "end time: " + opts.endT)
    else:
        opts.endT = None

    if opts.skip_clean and opts.skip_daily:
        parser.error(
            "Error: cannot skip both daily and clean averages")

    return (opts, indb)


def parse_localdata_for_comp(
        comp='Z', stdata=list, sta=None,
        start=UTCDateTime, end=UTCDateTime, ndval=nan):
    """
    Function to determine the path to data for a given component and alternate network

    Parameters
    ----------
    comp : str
        Channel for seismogram (one letter only)
    stdata : List
        Station list
    sta : Dict
        Station metadata from :mod:`~StDb` data base
    start : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time for request
    end : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time for request
    ndval : float or nan
        Default value for missing data

    Returns
    -------
    err : bool
        Boolean for error handling (`False` is associated with success)
    st : :class:`~obspy.core.Stream`
        Stream containing North, East and Vertical components of motion

    """

    from fnmatch import filter
    from obspy import read

    # Get start and end parameters
    styr = start.strftime("%Y")
    stjd = start.strftime("%j")
    edyr = end.strftime("%Y")
    edjd = end.strftime("%j")

    # Intialize to default positive error
    erd = True

    print(
        ("*          {0:2s}{1:1s} - Checking " +
            "Disk".format(sta.channel.upper(), comp.upper())))

    # Time Window Spans Single Day
    if stjd == edjd:
        # Format 1
        lclfiles = list(
            filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                   '.*.{4:2s}{5:1s}.SAC'.format(
                       styr, stjd, sta.network.upper(),
                       sta.station.upper(),
                       sta.channel.upper()[0:2],
                       comp.upper())))
        # Format 2
        if len(lclfiles) == 0:
            lclfiles = list(
                filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                       '.*.*{4:1s}.SAC'.format(
                           styr, stjd, sta.network.upper(),
                           sta.station.upper(), comp.upper())))
        # Alternate Nets (for CN/PO issues) Format 1
        if len(lclfiles) == 0:
            lclfiles = []
            for anet in sta.altnet:
                lclfiles.extend(
                    list(
                        filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                               '.*.{4:2s}{5:1s}.SAC'.format(
                                   styr, stjd, anet.upper(),
                                   sta.station.upper(),
                                   sta.channel.upper()[0:2],
                                   comp.upper()))))

        # Alternate Nets (for CN/PO issues) Format 2
        if len(lclfiles) == 0:
            # Check Alternate Networks
            lclfiles = []
            for anet in sta.altnet:
                lclfiles.extend(
                    list(
                        filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                               '.*.*{4:1s}.SAC'.format(
                                   styr, stjd, sta.network.upper(),
                                   sta.station.upper(), comp.upper()))))

        # If still no Local files stop
        if len(lclfiles) == 0:
            print("*              - Data Unavailable")
            return erd, None

        # Process the local Files
        for sacfile in lclfiles:
            # Read File
            st = read(sacfile, format="SAC")

            # Should only be one component, otherwise keep reading.
            # If more than 1 component, error
            if len(st) != 1:
                pass

            else:
                # Check for NoData and convert to NaN
                stnd = st[0].stats.sac['user9']
                eddt = False
                if (not stnd == 0.0) and (not stnd == -12345.0):
                    st[0].data[st[0].data == stnd] = ndval
                    eddt = True

                # Check start/end times in range
                if (st[0].stats.starttime <= start and
                        st[0].stats.endtime >= end):
                    st.trim(starttime=start, endtime=end)

                    # Check for Nan in stream
                    if True in isnan(st[0].data):
                        print(
                            "*          !!! Missing Data Present " +
                            "!!! Skipping (NaNs)")
                    else:
                        if eddt and (ndval == 0.0):
                            if any(st[0].data == 0.0):
                                print(
                                    "*          !!! Missing " +
                                    "Data Present !!! (Set to Zero)")

                        st[0].stats.update()
                        tloc = st[0].stats.location
                        if len(tloc) == 0:
                            tloc = "--"

                        # Processed succesfully...Finish
                        print(
                            ("*          {1:3s}.{2:2s}  - " +
                             "From Disk".format(
                                 st[0].stats.station,
                                 st[0].stats.channel.upper(), tloc)))
                        return False, st

    # Time Window spans Multiple days
    else:
        # Day 1 Format 1
        lclfiles1 = list(
            filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                   '.*.{4:2s}{5:1s}.SAC'.format(
                       styr, stjd, sta.network.upper(),
                       sta.station.upper(),
                       sta.channel.upper()[0:2],
                       comp.upper())))
        # Day 1 Format 2
        if len(lclfiles1) == 0:
            lclfiles1 = list(
                filter(
                    stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                    '.*.*{4:1s}.SAC'.format(
                        styr, stjd, sta.network.upper(),
                        sta.station.upper(), comp.upper())))
        # Day 1 Alternate Nets (for CN/PO issues) Format 1
        if len(lclfiles1) == 0:
            lclfiles1 = []
            for anet in sta.altnet:
                lclfiles1.extend(
                    list(
                        filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                               '.*.{4:2s}{5:1s}.SAC'.format(
                                   styr, stjd, anet.upper(),
                                   sta.station.upper(),
                                   sta.channel.upper()[0:2],
                                   comp.upper()))))
        # Day 1 Alternate Nets (for CN/PO issues) Format 2
        if len(lclfiles1) == 0:
            lclfiles1 = []
            for anet in sta.altnet:
                lclfiles1.extend(
                    list(
                        filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                               '.*.*{4:1s}.SAC'.format(
                                   styr, stjd, anet.upper(),
                                   sta.station.upper(),
                                   comp.upper()))))

        # Day 2 Format 1
        lclfiles2 = list(
            filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                   '.*.{4:2s}{5:1s}.SAC'.format(
                       edyr, edjd, sta.network.upper(
                       ), sta.station.upper(),
                       sta.channel.upper()[0:2],
                       comp.upper())))
        # Day 2 Format 2
        if len(lclfiles2) == 0:
            lclfiles2 = list(
                filter(stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                       '.*.*{4:1s}.SAC'.format(
                           edyr, edjd, sta.network.upper(),
                           sta.station.upper(), comp.upper())))
        # Day 2 Alternate Nets (for CN/PO issues) Format 1
        if len(lclfiles2) == 0:
            lclfiles2 = []
            for anet in sta.altnet:
                lclfiles2.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                            '.*.{4:2s}{5:1s}.SAC'.format(
                                edyr, edjd, anet.upper(),
                                sta.station.upper(),
                                sta.channel.upper()[0:2],
                                comp.upper()))))
        # Day 2 Alternate Nets (for CN/PO issues) Format 2
        if len(lclfiles2) == 0:
            lclfiles2 = []
            for anet in sta.altnet:
                lclfiles2.extend(
                    list(
                        filter(
                            stdata, '*/{0:4s}.{1:3s}.{2:s}.{3:s}' +
                            '.*.*{4:1s}.SAC'.format(
                                edyr, edjd, anet.upper(),
                                sta.station.upper(), comp.upper()))))

        # If still no Local files stop
        if len(lclfiles1) == 0 and len(lclfiles2) == 0:
            print("*              - Data Unavailable")
            return erd, None

        # Now try to merge the two separate day files
        if len(lclfiles1) > 0 and len(lclfiles2) > 0:
            # Loop over first day file options
            for sacf1 in lclfiles1:
                st1 = read(sacf1, format='SAC')
                # Loop over second day file options
                for sacf2 in lclfiles2:
                    st2 = read(sacf2, format='SAC')

                    # Check time overlap of the two files.
                    if (st1[0].stats.endtime >= st2[0].stats.starttime -
                            st2[0].stats.delta):
                        # Check for NoData and convert to NaN
                        st1nd = st1[0].stats.sac['user9']
                        st2nd = st2[0].stats.sac['user9']
                        eddt1 = False
                        eddt2 = False
                        if (not st1nd == 0.0) and (not st1nd == -12345.0):
                            st1[0].data[st1[0].data == st1nd] = ndval
                            eddt1 = True
                        if (not st2nd == 0.0) and (not st2nd == -12345.0):
                            st2[0].data[st2[0].data == st2nd] = ndval
                            eddt2 = True

                        st = st1 + st2
                        # Need to work on this HERE (AJS OCT 2015).
                        # If Calibration factors are different,
                        try:
                                # then the traces cannot be merged.
                            st.merge()

                            # Should only be one component, otherwise keep
                            # reading If more than 1 component, error
                            if len(st) != 1:
                                print(st)
                                print("merge failed?")

                            else:
                                if (st[0].stats.starttime <= start
                                        and st[0].stats.endtime >= end):
                                    st.trim(starttime=start, endtime=end)

                                    # Check for Nan in stream
                                    if True in isnan(st[0].data):
                                        print(
                                            "*          !!! Missing Data " +
                                            "Present !!! Skipping (NaNs)")
                                    else:
                                        if (eddt1 or eddt2) and (ndval == 0.0):
                                            if any(st[0].data == 0.0):
                                                print(
                                                    "*          !!! Missing " +
                                                    "Data Present !!! (Set " +
                                                    "to Zero)")

                                        st[0].stats.update()
                                        tloc = st[0].stats.location
                                        if len(tloc) == 0:
                                            tloc = "--"

                                        # Processed succesfully...Finish
                                        print(
                                            ("*          " +
                                             "{1:3s}.{2:2s}  - From " +
                                             "Disk".format(
                                                 st[0].stats.station,
                                                 st[0].stats.channel.upper(),
                                                 tloc)))
                                        return False, st

                        except:
                            pass
                    else:
                        print(("*                 - Merge Failed: " +
                               "No Overlap {0:s} - " +
                               "{1:s}".format(st1[0].stats.endtime,
                                              st2[0].stats.starttime -
                                              st2[0].stats.delta)))

    # If we got here, we did not get the data.
    print("*              - Data Unavailable")
    return erd, None


def get_data_NEZ(
        client=None, sta=None, start=UTCDateTime,
        end=UTCDateTime, stdata=list, ndval=nan):
    """
    Function to build a stream object for a seismogram in a given time window either
    by downloading data from the client object or alternatively first checking if the
    given data is already available locally.

    Note 
    ----
    Currently only supports NEZ Components!

    Parameters
    ----------
    client : :class:`~obspy.client.fdsn.Client`
        Client object
    sta : Dict
        Station metadata from :mod:`~StDb` data base
    start : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time for request
    end : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time for request
    stdata : List
        Station list
    ndval : float or nan
        Default value for missing data

    Returns
    -------
    err : bool
        Boolean for error handling (`False` is associated with success)
    trN : :class:`~obspy.core.Trace`
        Trace of North component of motion
    trE : :class:`~obspy.core.Trace`
        Trace of East component of motion
    trZ : :class:`~obspy.core.Trace` 
        Trace of Vertical component of motion

    """

    from fnmatch import filter
    from obspy import read
    from os.path import dirname, join, exists
    from numpy import any

    # Output
    print(("*     {0:s}.{1:2s} - NEZ:".format(sta.station,
                                              sta.channel.upper())))

    # Set Error Default to True
    erd = True

    # Check if there is local data
    if len(stdata) > 0:
        # Only a single day: Search for local data
        # Get Z localdata
        errZ, stZ = parse_localdata_for_comp(
            comp='Z', stdata=stdata, sta=sta, start=start, end=end,
            ndval=ndval)
        # Get N localdata
        errN, stN = parse_localdata_for_comp(
            comp='N', stdata=stdata, sta=sta, start=start, end=end,
            ndval=ndval)
        # Get E localdata
        errE, stE = parse_localdata_for_comp(
            comp='E', stdata=stdata, sta=sta, start=start, end=end,
            ndval=ndval)
        # Retreived Succesfully?
        erd = errZ or errN or errE
        if not erd:
            # Combine Data
            st = stZ + stN + stE

    # No local data? Request using client
    if erd:
        erd = False

        for loc in sta.location:
            tloc = loc
            # Constract location name
            if len(tloc) == 0:
                tloc = "--"
            # Construct Channel List
            channels = sta.channel.upper() + 'E,' + sta.channel.upper() + \
                'N,' + sta.channel.upper() + 'Z'
            print(("*          {1:2s}[ENZ].{2:2s} - Checking Network".format(
                sta.station, sta.channel.upper(), tloc)))
            try:
                st = client.get_waveforms(
                    network=sta.network,
                    station=sta.station, location=loc,
                    channel=channels, starttime=start,
                    endtime=end, attach_response=False)
                if len(st) == 3:
                    erd = False
                    print("*             - Data Downloaded")
                else:
                    erd = True
                    st = None
                    print("*             - Data Unavailable")
            except:
                erd = True
                st = None
                print("*             - Data Unavailable")

            # Break if we successfully obtained 3 components in st
            if not erd:
                break

    # Check the correct 3 components exist
    if st is None:
        print("* Error retrieving waveforms")
        print("****************************************************")
        return True, None, None, None
    elif (not st.select(component='Z')[0] or not st.select(component='E'[0])
          or not st.select(component='N')[0]):
        print("* Error retrieving waveforms")
        if not st.select(component='Z')[0]:
            print("*   --Z Missing")
        if not st.select(component='E')[0]:
            print("*   --E Missing")
        if not st.select(component='N')[0]:
            print("*   --N Missing")
        print("****************************************************")
        return True, None, None, None

    # Three components successfully retrieved
    print("* Waveforms Retrieved...")

    # Check trace lengths
    ll0 = len(st.select(component='Z')[0].data)
    ll1 = len(st.select(component='E')[0].data)
    ll2 = len(st.select(component='Z')[0].data)
    if not (ll0 == ll1 and ll0 == ll2):
        print(("* Error:  Trace lengths (Z,E,N): ", ll0, ll1, ll2))
        print("****************************************************")
        return True, None, None, None

    # Detrend data
    st.detrend('demean')
    st.detrend('linear')

    # Filter Traces
    st.filter('lowpass', freq=0.5*5.0, corners=2, zerophase=True)
    st.resample(5.0)
    st.filter('bandpass', freqmin=0.01, freqmax=0.25,
              corners=2, zerophase=True)

    # Extract traces
    trE = st.select(component='E')[0]
    trN = st.select(component='N')[0]
    trZ = st.select(component='Z')[0]

    # Return Flag and Data
    return False, trN, trE, trZ
