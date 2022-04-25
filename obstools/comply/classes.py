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


import numpy as np
import pickle
from obstools.atacr import utils, DayNoise, StaNoise
np.seterr(all='ignore')


class Comply(object):
    """
    A Comply object contains attributes that calculate and store the
    compliance and coherence functions for the available channels.

    Note
    ----
    The object is initialized with either a processed
    :class:`~obstools.atacr.classes.DayNoise` or
    :class:`~obstools.atacr.classes.StaNoise` object. Each individual
    spectral quantity is unpacked as an object attribute, but all of them
    are discarded as the object is saved to disk.

    Attributes
    ----------
    elev : float
        Station elevation in meters (OBS stations have negative elevations)
    f : :class:`~numpy.ndarray`
        Frequency axis for corresponding time sampling parameters
    c11 : `numpy.ndarray`
        Power spectra for component `H1`. Other identical attributes are
        available for the power, cross and rotated spectra:
        [11, 12, 1Z, 1P, 22, 2Z, 2P, ZZ, ZP, PP, HH, HZ, HP]
    tf_list : Dict
        Dictionary of possible transfer functions from the available
        components (obtained from the
        :class:`~obstools.atacr.classes.DayNoise` or the
        :class:`~obstools.atacr.classes.StaNoise` noise objects)
    complyfunc : Dict
        Dictionary of compliance and coherence functions given the available
        components.

    """

    def __init__(self, objnoise=None, elev=None):

        if any(value is None for value in [elev, objnoise]):
            raise(Exception(
                "Error: Initializing EventStream object with None values - " +
                "aborting"))

        if (objnoise and not isinstance(objnoise, DayNoise) and
                not isinstance(objnoise, StaNoise)):
            msg = "Error: A TFNoise object must be initialized with only " +\
                "one of type DayNoise or StaNoise object"
            raise TypeError(msg)

        if not objnoise.av:
            raise(Exception(
                "Error: Noise object has not been processed (QC and " +
                "averaging) - aborting"))

        self.elevation = elev
        self.f = objnoise.f
        self.c11 = objnoise.power.c11
        self.c22 = objnoise.power.c22
        self.cZZ = objnoise.power.cZZ
        self.cPP = objnoise.power.cPP
        self.cHH = objnoise.rotation.cHH
        self.cHZ = objnoise.rotation.cHZ
        self.cHP = objnoise.rotation.cHP
        self.c12 = objnoise.cross.c12
        self.c1Z = objnoise.cross.c1Z
        self.c1P = objnoise.cross.c1P
        self.c2Z = objnoise.cross.c2Z
        self.c2P = objnoise.cross.c2P
        self.cZP = objnoise.cross.cZP
        self.tf_list = objnoise.tf_list

    class ComplyDict(dict):

        def __init__(self):
            self = dict()

        def add(self, key, value):
            self[key] = value

    def calculate_compliance(self):
        """
        Method to calculate compliance and coherence functions from the
        averaged (daily or station-averaged) noise spectra.

        Attributes
        ----------
        complyfunc : Dict
            Container Dictionary for all possible compliance and
            coherence functions

        Examples
        --------

        Calculate compliance and coherence functions for a DayNoise object.
        In these examples, station elevation is extracted from the IRIS
        metadata aggregator site: http://ds.iris.edu/mda/7D/M08A/

        >>> from obstools.atacr import DayNoise
        >>> from obstools.comply import Comply
        >>> daynoise = DayNoise('demo')
        Uploading demo data - March 04, 2012, station 7D.M08A
        >>> daynoise.QC_daily_spectra()
        >>> daynoise.average_daily_spectra()
        >>> daycomply = Comply(objnoise=daynoise, elev=-126.4)
        >>> daycomply.calculate_compliance()
        >>> tfnoise.complyfunc.keys()
        dict_keys(['ZP', 'ZP-21', 'ZP-H'])

        Calculate compliance and coherence functions for a StaNoise object

        >>> from obstools.atacr import StaNoise
        >>> from obstools.comply import Comply
        >>> stanoise = StaNoise('demo')
        Uploading demo data - March 01 to 04, 2012, station 7D.M08A
        >>> stanoise.QC_sta_spectra()
        >>> stanoise.average_sta_spectra()
        >>> stacomply = Comply(objnoise=stanoise, elev=-126.4)
        >>> stacomply.calculate_compliance()
        >>> stacomply.complyfunc.keys()
        dict_keys(['ZP', 'ZP-21'])

        """

        def wavenumber(omega, H):
            """
            Function to approximate wavenumber from dispersion relation

            H is depth below the seafloor, in meters
            omega is a vector of positive angular frequencies

            Stephen G. Mosher, 2020

            """

            import numpy.polynomial as poly

            g = 9.79329
            N = len(omega)

            # Approximations for k when k*H is very large (deep case) or
            # very small (shallow case)
            k_deep = omega**2 / g
            k_shal = omega / np.sqrt(g * H)

            """
            Alternatively, we can use a rational approximation to
            tanh(x) to solve k for any omega. This approximation gives
            a quartic equation, we take the positive real roots as the
            value of k we're interested in. The rational approximation
            being used is always better than the shallow approximation.
            However, it's only better than the deep approximation if
            k*H < 2.96. Therefore, we keep the solutions to k we find,
            using the rational approximation for k*H < 2.96 and use the
            deep water approximation to solve for k otherwise. The
            average error is just under 1% and the maximum error is
            2.5%.
            """

            k = np.zeros(len(omega))

            for i, om in enumerate(omega):

                if i == 0:
                    k[i] = 0.
                else:

                    a0 = -27 * om**2 / g    # constant terms
                    a1 = 0.                 # no linear terms
                    a2 = 27 * H - (9 * om**2 * H**2)/g     # quadratic terms
                    a3 = 0.                 # no cubic terms
                    a4 = H**3           # quartic terms

                    p = poly.Polynomial([a0, a1, a2, a3, a4])
                    solu = poly.Polynomial.roots(p)
                    positive_roots = solu[solu > 0]
                    real_positive_root = \
                        positive_roots[positive_roots.imag == 0].real[0]
                    k[i] = real_positive_root

            # For k*H >= 2.96, prefer the deep approximation above
            for i, wavenumber in enumerate(k_deep):
                if wavenumber * H > 2.96:
                    k[i] = k_deep[i]

            return k

        # Calculate wavenumber - careful here, elevation is negative
        k = wavenumber(2.*np.pi*self.f, -1.*self.elevation)

        # Initialize empty dictionary
        complyfunc = self.ComplyDict()

        # Cycle through all available transfer functions in the objnoise
        # object
        for key, value in self.tf_list.items():

            if key == 'ZP':
                if value:
                    admit_ZP = utils.admittance(self.cZP, self.cPP)
                    compl_ZP = k*admit_ZP
                    coh_ZP = utils.coherence(self.cZP, self.cPP, self.cZZ)
                    complyfunc.add('ZP', [compl_ZP, coh_ZP])

            elif key == 'ZP-21':
                if value:
                    lc1cZ = np.conj(self.c1Z)/self.c11
                    lc1c2 = np.conj(self.c12)/self.c11
                    lc1cP = np.conj(self.c1P)/self.c11

                    coh_12 = utils.coherence(self.c12, self.c11, self.c22)
                    coh_1P = utils.coherence(self.c1P, self.c11, self.cPP)
                    coh_1Z = utils.coherence(self.c1Z, self.c11, self.cZZ)

                    gc2c2_c1 = self.c22*(1. - coh_12)
                    gcPcP_c1 = self.cPP*(1. - coh_1P)
                    gcZcZ_c1 = self.cZZ*(1. - coh_1Z)

                    gc2cZ_c1 = np.conj(self.c2Z) - np.conj(lc1c2*self.c1Z)
                    gcPcZ_c1 = self.cZP - np.conj(lc1cP*self.c1Z)

                    gc2cP_c1 = np.conj(self.c2P) - np.conj(lc1c2*self.c1P)

                    lc2cP_c1 = gc2cP_c1/gc2c2_c1
                    lc2cZ_c1 = gc2cZ_c1/gc2c2_c1

                    coh_c2cP_c1 = utils.coherence(gc2cP_c1, gc2c2_c1,
                                                  gcPcP_c1)
                    coh_c2cZ_c1 = utils.coherence(gc2cZ_c1, gc2c2_c1,
                                                  gcZcZ_c1)

                    gcPcP_c1c2 = gcPcP_c1*(1. - coh_c2cP_c1)
                    gcPcZ_c1c2 = gcPcZ_c1 - np.conj(lc2cP_c1)*gc2cZ_c1
                    gcZcZ_c1c2 = gcZcZ_c1*(1. - coh_c2cZ_c1)

                    admit_ZP_21 = utils.admittance(
                        gcPcZ_c1c2, gcPcP_c1c2)
                    compl_ZP_21 = k*admit_ZP_21
                    coh_ZP_21 = utils.coherence(
                        gcPcZ_c1c2, gcPcP_c1c2, gcZcZ_c1c2)

                    complyfunc.add('ZP-21', [compl_ZP_21, coh_ZP_21])

            elif key == 'ZP-H':
                if value:
                    lcHcP = np.conj(self.cHP)/self.cHH

                    coh_HP = utils.coherence(self.cHP, self.cHH, self.cPP)
                    coh_HZ = utils.coherence(self.cHZ, self.cHH, self.cZZ)

                    gcPcP_cH = self.cPP*(1. - coh_HP)
                    gcZcZ_cH = self.cZZ*(1. - coh_HZ)
                    gcPcZ_cH = self.cZP - np.conj(lcHcP*self.cHZ)

                    admit_ZP_H = utils.admittance(gcPcZ_cH, gcPcP_cH)
                    compl_ZP_H = k*admit_ZP_H
                    coh_ZP_H = utils.coherence(gcPcZ_cH, gcPcP_cH, gcZcZ_cH)

                    complyfunc.add('ZP-H', [compl_ZP_H, coh_ZP_H])

            self.complyfunc = complyfunc

    def save(self, filename, form='pkl'):
        """
        Method to save the object to file using `~Pickle`.

        Parameters
        ----------
        filename : str
            File name

        Examples
        --------

        Following from the example outlined in method
        :func:`~obstools.comply.classes.Comply.calculate_compliance`, we simply
        save the final object

        >>> daycomply.save('daycomply_demo.pkl')

        Check that object has been saved

        >>> import glob
        >>> glob.glob("./daycomply_demo.pkl")
        ['./daycomply_demo.pkl']

        """

        if not self.complyfunc:
            print("Warning: saving before having calculated the compliance " +
                  "and coherence functions")

        if form == 'pkl':

            # Remove traces to save disk space
            del self.c11
            del self.c22
            del self.cZZ
            del self.cPP
            del self.cHH
            del self.cHZ
            del self.cHP
            del self.c12
            del self.c1Z
            del self.c1P
            del self.c2Z
            del self.c2P
            del self.cZP

            file = open(filename.parent / (filename.name + '.' + form), 'wb')
            pickle.dump(self, file)
            file.close()

        elif form == 'csv':

            import pandas as pd

            if 'ZP-H' in self.complyfunc:
                df = pd.DataFrame(
                    {'Frequency': self.f,
                     'Compliance ZP': self.complyfunc['ZP'][0],
                     'Coherence ZP': self.complyfunc['ZP'][1],
                     'Compliance ZP-21': self.complyfunc['ZP-21'][0],
                     'Coherence ZP-21': self.complyfunc['ZP-21'][1],
                     'Compliance ZP-H': self.complyfunc['ZP-H'][0],
                     'Coherence ZP-H': self.complyfunc['ZP-H'][1]})
            elif 'ZP-21' in self.complyfunc:
                df = pd.DataFrame(
                    {'Frequency': self.f,
                     'Compliance ZP': self.complyfunc['ZP'][0],
                     'Coherence ZP': self.complyfunc['ZP'][1],
                     'Compliance ZP-21': self.complyfunc['ZP-21'][0],
                     'Coherence ZP-21': self.complyfunc['ZP-21'][1]})
            else:
                df = pd.DataFrame(
                    {'Frequency': self.f,
                     'Compliance ZP': self.complyfunc['ZP'][0],
                     'Coherence ZP': self.complyfunc['ZP'][1]})

            df.to_csv(filename.parent / (filename.name +
                      '.' + form), index=False)

        return
