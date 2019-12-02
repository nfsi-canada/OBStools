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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

ATaCR is a module for the correction of vertical component data from OBS
stations from tilt and compliance noise. This module is a translation of the 
Matlab code `ATaCR  <https://github.com/helenjanisz/ATaCR>`_ and the acronym 
stands for Automatic Tilt and Compliance Removal. 

The analysis can be carried out for either one (or both) compliance or tilt
corrections. In all cases the analysis requires at least vertical component
data. Additional data required depend on the type of analysis. The software 
will automatically calculate all possible corrections depending on the available
channels. For more details on the 
theory and methodology, we refer the interested reader to the following papers:

- Bell, S. W., D. W. Forsyth, and Y. Ruan (2014), Removing noise from the 
  vertical component records of ocean-bottom seismometers: Results from year one of the 
  Cascadia Initiative, Bull. Seismol. Soc. Am., 105, 300-313.

- Crawford, W.C., Webb, S.C., (2000). Identifying and removing tilt noise from 
  low-frequency (0.1 Hz) seafloor vertical seismic data, Bull. seism. Soc. Am., 90, 952-963.

- Janiszewski, H A, J B Gaherty, G A Abers, H Gao, Z C Eilon, Amphibious surface-wave 
  phase-velocity measurements of the Cascadia subduction zone, Geophysical Journal 
  International, Volume 217, Issue 3, June 2019, Pages 1929-1948, https://doi.org/10.1093/gji/ggz051



Corrections
***********

**Compliance**

Compliance is defined as the spectral ratio between pressure and vertical
displacement data. Compliance noise arises from seafloor deformation due
to seafloor and water wave effects including infragravity waves). 
This is likely the main source of noise in vertical component OBS data. 
This analysis therefore requires both vertical (?HZ) and pressure (?XH) data.

**Tilt**

Tilt noise arises from OBS stations that are not perfectly leveled, and
therefore the horizontal seafloor deformation leaks into the vertical
component. This effect can be removed by calculating the spectral
ratio between horizontal and vertical displacement data. In most cases,
however, the tilt direction (measured on a compass - as opposed to tilt
angle, measured from the vertical axis) is unknown and must be determined
from the coherence between rotated horizontal components and the vertical
component. This analysis therefore requires vertical (?HZ) and the two
horizontal (?H1,2) component data.

**Compliance + Tilt**

It is of course possible to combine both corrections and apply them
sequentially. In this case the tilt noise is removed first, then compliance.
This analysis requires all four components: three-component
seismic (?HZ,1,2) and pressure (?XH) data.

"""