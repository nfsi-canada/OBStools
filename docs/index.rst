.. figure:: ../stdb/examples/figures/StDb_logo.png
   :align: center

OBStools - documentation
========================

OBStools is a package containing Python tools for processing broadband
ocean-bottom seismic (OBS) data. Current functionalities include removing
vertical component noise from tilt and compliant effects. The code uses 
the ``StDb`` package for querying and building a station database and
is used through command-line scripts.

The software is a Python translation of the Matlab software 
`ATaCR <https://github.com/helenjanisz/ATaCR>`_.


Quick links
"""""""""""

* `OBStools Git repository <https://github.com/paudetseis/OBStools>`_
* `StDb Git repository <https://github.com/schaefferaj/StDb>`_
* `ATaCR Git repository <https://github.com/helenjanisz/ATaCR>`_

.. Getting Started
.. """""""""""""""

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   api
   database

.. Content
.. """""""

.. toctree::
   :maxdepth: 1
   :caption: Content

   classes
   io
   convert
   gui
   kml

.. Programs & Tutorials
.. """"""""""""""""""""

.. toctree::
   :maxdepth: 1
   :caption: Programs & Tutorials

   query
   ls
   merge
   gen
   append
   edit
   dump
   stdb2kml
