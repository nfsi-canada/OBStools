# Copyright 2019 Pascal Audet
#
# This file is part of RfPy.
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

Module containing the main utility functions used in the `RfPy` scripts
that accompany this package.

"""

# -*- coding: utf-8 -*-
from obspy import UTCDateTime
from numpy import nan, isnan
from obspy.core import Stream, read


def get_orient_options():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    from optparse import OptionParser, OptionGroup
    from os.path import exists as exist
    from obspy import UTCDateTime
    from numpy import nan

    parser = OptionParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to find orientation of station from "+
        "receiver function data ")

    # General Settings
    parser.add_option(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys for " +
        "which to perform the analysis. These must be " +
        "contained within the station database. Partial keys will " +
        "be used to match against those in the dictionary. For " +
        "instance, providing IU will match with all stations in " +
        "the IU network [Default processes all stations in the database]")
    parser.add_option(
        "-v", "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")
    parser.add_option(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing figures. " +
        "[Default False]")

    PreGroup = OptionGroup(
        parser,
        title='Pre-processing Settings',
        description="Options for pre-processing of receiver function " +
        "data before plotting")
    PreGroup.add_option(
        "--snr",
        action="store",
        type=float,
        dest="snr",
        default=-9999.,
        help="Specify the SNR threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_option(
        "--snrh",
        action="store",
        type=float,
        dest="snrh",
        default=-9999.,
        help="Specify the horizontal component SNR threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_option(
        "--cc",
        action="store",
        type=float,
        dest="cc",
        default=-1.,
        help="Specify the CC threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_option(
        "--no-outlier",
        action="store_true",
        dest="no_outl",
        default=False,
        help="Set this option to delete outliers based on the MAD "+
        "on the variance. [Default False]")
    PreGroup.add_option(
        "--bp",
        action="store",
        type=str,
        dest="bp",
        default=None,
        help="Specify the corner frequencies for the bandpass filter. " +
        "[Default no filtering]")
    PreGroup.add_option(
        "--pws",
        action="store_true",
        dest="pws",
        default=False,
        help="Set this option to use phase-weighted stacking during binning "+
        " [Default False]")
    PreGroup.add_option(
        "--nbaz",
        action="store",
        dest="nbaz",
        type=int,
        default=72,
        help="Specify integer number of back-azimuth bins to consider " +
        "(typically 36 or 72). If not None, the plot will show receiver " +
        "functions sorted by back-azimuth values. [Default 72]")
    PreGroup.add_option(
        "--trange",
        action="store",
        default=None,
        type=str,
        dest="trange",
        help="Specify the time range for decomposition (sec). Negative times "+
        "are allowed [Default -1., 1.]")
    PreGroup.add_option(
        "--boot",
        action="store_true",
        dest="boot",
        default=False,
        help="Set this option to calculate bootstrap statistics "+
        " [Default False]")
    PreGroup.add_option(
        "--plot-f",
        action="store_true",
        dest="plot_f",
        default=False,
        help="Set this option to plot the function f(phi) "+
        "[Default False]")
    PreGroup.add_option(
        "--plot-comps",
        action="store_true",
        dest="plot_comps",
        default=False,
        help="Set this option to plot the misoriented and rotated harmonic "+
        "components [Default False]")




    parser.add_option_group(PreGroup)

    (opts, args) = parser.parse_args()

    # Check inputs
    if len(args) != 1:
        parser.error("Need station database file")
    indb = args[0]
    if not exist(indb):
        parser.error("Input file " + indb + " does not exist")

    # create station key list
    if len(opts.stkeys) > 0:
        opts.stkeys = opts.stkeys.split(',')

    if opts.bp is not None:
        opts.bp = [float(val) for val in opts.bp.split(',')]
        opts.bp = sorted(opts.bp)
        if (len(opts.bp)) != 2:
            parser.error(
                "Error: --bp should contain 2 " +
                "comma-separated floats")

    if opts.trange is None:
        opts.tmin = -1.
        opts.tmax = 1.
    if opts.trange is not None:
        opts.trange = [float(val) for val in opts.trange.split(',')]
        opts.trange = sorted(opts.trange)
        if (len(opts.trange)) != 2:
            parser.error(
                "Error: --trange should contain 2 " +
                "comma-separated floats")

    return (opts, indb)

