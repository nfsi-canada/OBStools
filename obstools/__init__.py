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

Licence
-------

Copyright 2019 Pascal Audet

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Installation
------------

Dependencies
++++++++++++

- `StDb <https://github.com/schaefferaj/StDb>`_
- `ObsPy <https://github.com/obspy/obspy>`_

Conda environment
+++++++++++++++++

We recommend creating a custom ``conda`` environment
where ``OBStools`` can be installed along with its dependencies.

.. sourcecode:: bash

   conda create -n obs -c conda-forge python=3.12 obspy

Activate the newly created environment:

.. sourcecode:: bash

   conda activate obs

Install StDb:

.. sourcecode:: bash

   pip install git+https://github.com/schaefferaj/stdb

Installing development branch on GitHub
+++++++++++++++++++++++++++++++++++++++

.. sourcecode:: bash

   pip install git+https://github.com/nfsi-canada/obstools

Installing from source
++++++++++++++++++++++

- Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/paudetseis/OBStools.git
   cd OBStools

- Install using pip:

.. sourcecode:: bash

   pip install .

Using local data
----------------

The various scripts packaged with ``OrientPy`` use FDSN web services
through and ``ObsPy`` `Client` to load waveform data. For waveform
data locally stored on your hard drive, the scripts can use a `Client`
that reads a `SeisComP Data Structure 
<https://docs.obspy.org/packages/autogen/obspy.clients.filesystem.sds.html>`_ 
archive containing SAC or miniSEED waveform data. Check out the scripts 
``bng_calc`` and ``dl_calc`` below and the argument ``--SDS-path`` and 
``--dtype`` for more details.

Station Metadata
++++++++++++++++

If you have data stored locally on your drive, it is likely you also
have a station `XML <https://www.fdsn.org/xml/station/>`_ file
containing the metadata. The corresponding ObsPy documentation is
`here <https://docs.obspy.org/packages/obspy.core.inventory.html>`_. 

You can now use a stationXML (`.xml`) file instead of the StDb `.pkl` format. 
Alternatively, you can convert the stationXML file to an StDb `.pkl` file
by running the command ``gen_stdb station.xml`` (these options are only
available on StDb version 0.2.7. If you don't have a station `XML` file but you have
a dataless SEED file, you can convert it first to XML using `this
tools <https://seiscode.iris.washington.edu/projects/stationxml-converter>`_.

Waveform Data
+++++++++++++

The SDS folder containing the waveform data has the structure:

.. code-block:: python

   archive
     + year
       + network code
         + station code
           + channel code + type
             + one file per day and location, e.g. NET.STA.LOC.CHAN.TYPE.YEAR.DOY


For example:

.. code-block:: python

   SDS/
     2014/
       YH/
         LOBS3/
           HH1.D/ 
             YH.LOBS3..CH1.D.2014.332
             ...


Note, the filename does not include the extension (`.MSEED` or `.SAC`),
and the characters `.D` (for type Data) that appear in both the
channel code and the filename. Note also the two dots (`..`). If
there is a location code, it should appear between those dots (e.g.,
for a location code `10`, the corresponding filename should be
`YH.LOBS3.10.HH1.D.2014.332`). There is no location code for the 
YH.LOBS3 data, and this field is simply absent from the filenames.
Finally, the day-of-year (DOY) field must be zero-padded to be exactly
3 characters.
"""

__version__ = '0.2.1'

__author__ = 'Pascal Audet'

from . import atacr
from . import comply
