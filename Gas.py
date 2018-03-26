#Embedded file name: /scratch/users/alatimer/meoh-vs-methane/analysis/catclass.py
from ase import Atoms
from ase.io import read, write
import pickle
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from glob import glob
import re

class Gas:
    def __init__(self, 
            name = None, 
            symmetrynumber = None,
            geometry = None,
            spin = None,
            ):

        self.name = name
        self.symmetrynumber = symmetrynumber
        self.geometry = geometry
        self.spin = spin


class GasCalc:
    """
    """

    def __init__(self, 
            name = None, 
            psp_dict = None,
            pw = None,
            dw = None,
            xc = None,
            energy = None,
            ):

        self.name = name
        self.psp_dict = psp_dict
        self.pw = pw
        if dw == None:
            self.dw = pw*10
        self.xc = xc
        self.energy = energy

class GasCalcs:
    """
    """

    def __init__(self, classes):
        self.data = classes

    def filter(self, fun, *args, **kwargs):
        """
        Takes list of chargexfer objects and filters out those that do not give fun(chargexfer) == True.
        Returns new ChargeXfers object
        """
        out = []
        for c in self.data:
            if fun.__name__ == '<lambda>':
                bool = fun(c)
            else:
                bool = fun(c, *args, **kwargs)
            if bool:
                out.append(c)

        return classes(out)

    def get_property(self, ppt):
        """
        Return list of property for each chargexfer
        """
        out = []
        for c in self.data:
            eval('out.append(c.%s)' % ppt)
        return np.array(out)
