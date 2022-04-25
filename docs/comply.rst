Comply
======

Introduction
************

``comply`` is a module for the calculation of compliance and coherence 
functions from OBS stations. These functions can be used to retrieve shear-wave
velocity structure in the shallow subsurface. The compliance and coherence functions
use either :class:`~obstools.atacr.classes.DayNoise` or :class:`~obstools.atacr.classes.StaNoise`
objects as input. If horizontal components are available in these objects, tilt
is removed prior to compliance and coherence calculations. For more details on the 
theory and applications, we refer the interested reader to the following papers:

- Crawford, W.C., S. Webb, S.C, and Hildebrand, J. (1991). Sea-floor compliance 
  observed by long-period pressure and displacement measurements, J. Geophys. Res., 
  96(B10), 16,151–16,160.

- Doran, A., and Laske, G. (2019). Seismic structure of marine sediments and upper 
  oceanic crust surrounding Hawaii, J. Geophys. Res., 124, 2038-2056, 
  https://doi.org/10.1029/2018JB016548.

- Zha, Y., and Webb, S.C. (2016). Crustal shear velocity structure in the Southern Lau
  Basin constrained by seafloor compliance, J. Geophys.  Res., 121, 3220-3237,
  https://doi.org/10.1002/2015JB012688.

Compliance
**********

Compliance is defined as the spectral ratio between pressure and vertical
displacement data. Compliance arises from vertical seafloor deformation due
to infragravity wave effects, which propagate through the deep 
oceans with periods longer than 30 seconds. The normalized compliance
:math:`\eta(\omega)` is defined as

.. math::

   \eta(\omega) = k(\omega)\frac{u_z(\omega)}{p(\omega)}

where :math:`\omega` is the angular frequency, :math:`k` is the wavenumber,
:math:`u_z` is the vertical displacement spectrum and :math:`p` is the 
pressure spectrum. The wavenumber :math:`k` is obtained by solving the 
dispersion relation

.. math::

   \omega = gk(\omega)\tanh(k(\omega)H)

where :math:`g` is the gravitational acceleration and :math:`H` is
the seafloor depth (positive downward). 

Tilting of the OBS sensor due to seafloor currents can obfuscate the
coherence between vertical and pressure components in the infra-gravity band.
Tilt noise is therefore removed from the vertical component data prior
to calculating compliance.

API documentation
*****************

:mod:`~obstools.comply` defines the following class:

- :class:`~obstools.comply.classes.Comply`

The class :class:`~obstools.comply.classes.Comply` contains attributes
and methods for the calculation of the compliance and coherence functions 
from noise traces. A ``Comply`` object works with 
either one of :class:`~obstools.atacr.classes.DayNoise` and 
:class:`~obstools.atacr.classes.StaNoise` objects to calculate all possible
compliance and coherence functions across all available components. 
These transfer functions are saved as attributes of the object in a Dictionary. 

.. note::

    In the examples below, the SAC data were obtained and pre-processed
    using the accompanying scripts ``atacr_download_data``. See the 
    script and tutorial for details.

.. autoclass:: obstools.comply.classes.Comply
   :members:

Scripts
*******

There is only one Python scripts that accompanies ``~obstools.comply``. This script can be used
in bash scripts to automate data processing. Please see the :mod:`~obstools.atacr` 
module and scripts to prepare the data for use in the following script.

``comply_calculate``
++++++++++++++++++++++++++++++++++

Description
-----------

Calculates compliance and coherence functions using the noise windows flagged as *good*, for either
a single day (from ``atacr_daily_spectra``) or for those averaged over several days
(from ``atacr_clean_spectra``), if available. The transfer functions are stored to disk.
Station selection is specified by a network and station code. The database is 
provided as a ``StDb`` dictionary.

Usage
-----

.. code-block::

    $ comply_calculate -h
    usage: comply_calculate [options] <Station Database>

    Script used to calculate compliance functions between various components. The noise data can be those obtained
    from the daily spectra (i.e., from `atacr_daily_spectra`) or those obtained from the averaged noise spectra
    (i.e., from `atacr_clean_spectra`). Flags are available to specify the source of data to use as well as the
    time range over which to calculate the transfer functions. The stations are processed one by one and the data
    are stored to disk.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      --keys STKEYS         Specify a comma separated list of station keys for which to perform the analysis.
                            These must be contained within the station database. Partial keys will be used to
                            match against those in the dictionary. For instance, providing IU will match with all
                            stations in the IU network. [Default processes all stations in the database]
      -O, --overwrite       Force the overwriting of pre-existing data. [Default False]
      --save-format SAVEFORMAT
                            Specify the format of the output files. Options are: 'pkl' or 'csv'. [Default 'pkl']

    Time Search Settings:
      Time settings associated with searching for day-long seismograms

      --start STARTT        Specify a UTCDateTime compatible string representing the start day for the data
                            search. This will override any station start times. [Default start date of each
                            station in database]
      --end ENDT            Specify a UTCDateTime compatible string representing the start time for the data
                            search. This will override any station end times. [Default end date of each station in
                            database]

    Parameter Settings:
      Miscellaneous default values and settings

      --skip-daily          Skip daily spectral averages in construction of compliance and coherence functions.
                            [Default False]
      --skip-clean          Skip cleaned spectral averages in construction of compliance and coherence functions.
                            Defaults to True if data cannot be found in default directory. [Default False]

    Figure Settings:
      Flags for plotting figures

      --fig                 Plot compliance and coherence functions figure. [Default does not plot figure]
      --save-fig            Set this option if you wish to save the figure(s). [Default does not save figure]
      --format FORM         Specify format of figure. Can be any one of the validmatplotlib formats: 'png', 'jpg',
                            'eps', 'pdf'. [Default 'png']

Tutorial
********

.. note::

    Here we build on the steps already carried out in the ``ATaCR`` tutorial. Specifically,
    steps ``0`` through ``3`` should be done prior to performing the following step.

Compliance calculation
++++++++++++++++++++++

Once the ``StaNoise`` objects have been produced and saved to disk, the compliance and coherence 
functions across all available components can be calculated. By default the software
will calculate the ones for which the spectral averages are available. 

If no horizontal components are available (i.e., only ``?HZ`` and ``?DH?`` components are available),
the only compliance function possible is:

- ``ZP``

If horizontal components are available, the 
following compliance functions are possible:

- ``ZP``

- ``ZP-21``

If you are using a ``DayNoise`` object to calculate the compliance functions,
the following may also be possible (if all components are available):

- ``ZP-H``

In this example we calculate all available compliance functions for all available data.
In this case we do not need to specify any option and type in a terminal:

.. code-block::

    $ comply_calculate M08A.pkl

    Path to COMPL_STA/7D.M08A doesn't exist - creating it
     
     
    |===============================================|
    |===============================================|
    |                       M08A                    |
    |===============================================|
    |===============================================|
    |  Station: 7D.M08A                             |
    |      Channel: BH; Locations: --               |
    |      Lon: -124.90; Lat:  44.12                |
    |      Start time: 2011-10-20 00:00:00          |
    |      End time:   2012-07-18 23:59:59          |
    |-----------------------------------------------|

    ************************************************************
    * Calculating compliance functions for key 7D.M08A and day 2012.069

    ************************************************************
    * Calculating compliance functions for key 7D.M08A and day 2012.088

        ...

    **********************************************************************
    * Calculating compliance functions for key 7D.M08A and range 2011.293-2012.200.

Note how the ``DayNoise`` objects are read randomly from disk, followed by the 
``StaNoise`` object. The result is a ``Comply`` object that is saved to a newly
created folder called ``COMP_STA/7D.M08A/``. Once the objects are saved to file,
you should  use the option ``-O`` to overwrite any file on disk.

We can produce a Figure by re-running the previous command with the options
``--fig``. We also specify the low-frequency limit to be 0.02 Hz with the option
``--f0=0.02``, which we select based on the observed
non-physical kink in the compliance curves (you can see it by ignoring the
argument ``--f0``)

.. figure:: ../obstools/examples/figures/Figure_comply.png
   :align: center

Figure 1: Coherence and compliance functions for the component combinations 
of interest, as indicated in the title of each subplot. This example is for 
the month of March 2012 for station M08A. The QC'ed daily compliance functions 
are shown in grey. The QC'ed stations average calculated for the whole month is 
shown in black. Note that *all* daily compliance functions are calculated and 
displayed, even those that did not pass the QC threshold for the station average.

