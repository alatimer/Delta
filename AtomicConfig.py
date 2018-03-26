#Embedded file name: /scratch/users/alatimer/meoh-vs-methane/analysis/catclass.py
from ase import Atoms
from ase.io import read, write
import pickle
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from glob import glob
import re

class AtomicConfig:
    """
    """

    def __init__(self, 
            surf_name = None, 
            surf_class =  None,
            species_name = None, 
		    traj_loc = None, 
            vib_loc = None,
            beef_loc= None,
            species_type = None,
            symmetrynumber = None,
            geometry = None,
            spin = None,
            calc_params = None,
            tag = '',
            ):

        self.surf_name = surf_name
        self.surf_class = surf_class
        self.species_name = species_name
        self.species_type = species_type
        self.traj_loc = traj_loc
        self.symmetrynumber = symmetrynumber
        self.geometry = geometry
        self.spin = spin
        self.tag = tag

        if calc_params != None:
            self.calc_params = calc_params
        if traj_loc != None:
            self.atoms = read(traj_loc)
            self.calc_params = self.get_calc_params()
        if beef_loc != None:   
            self.beef = self.beef_reader()
        if vib_loc != None:
            self.vibs = self.vibs_reader(vib_loc)
            if self.species_type == 'gas':
                self.gibbs = IdealGasThermo(vib_energies=self.vibs,
                                potentialenergy=0,
                                atoms=self.atoms,
                                geometry=self.geometry,
                                symmetrynumber=self.symmetrynumber,
                                spin=self.spin,
                                )
            elif self.species_type == 'adsorbate':
                self.gibbs = HarmonicThermo(vib_energies = self.vibs, potentialenergy = 0)
            else:
                self.gibbs = None

    def get_calc_params(self):
        directory = '/'.join(self.traj_loc.split('/')[0:-1])
        if len(directory) < 1:
            directory = '.'
        calcdirs = glob(directory+'/ekspressen*/log')+glob(directory+'/calcdir/log')
        log_file = open(calcdirs[0],'r')
        inpdirs = glob(directory+'/ekspressen*/pw.inp')+glob(directory+'/calcdir/pw.inp')
        inp_file = open(inpdirs[0],'r')
        calc_params = {}
        calc_params['psp'] = {}
        for line in log_file:
            line = line.strip()
            #XC
            if line.startswith('Exchange-correlation'):
                if 'BEEF' in line:
                    xc='beef'
                elif 'RPBE' in line:
                    xc = 'rpbe'
                elif 'PBE' in line:
                    xc = 'pbe'
            #PSP
            if  line.startswith('PseudoPot.'):
                elem = line.split()[4]
            if line.startswith('MD5 check sum:'):
                md5 = line.split()[-1]
                calc_params['psp'][elem]=md5
            #PW
            if line.startswith('kinetic-energy cutoff'):
                pw = float(re.findall("\d+\.\d+", line)[0])
                pw *=13.61 #ryd to ev
                pw = str(int(pw)) #round pw and turn to str
        calc_params['xc'] = xc
        calc_params['pw'] = pw
        #Kpts
        inplines = inp_file.readlines()
        kpts = inplines[-1].split()
        calc_params['kpts'] = ''.join(kpts[0:3])
        
        if len(calc_params)<3:
            print "Unable to extract calculator parameters automatically, user can supply manually."
            exit()
        return calc_params

    def beef_reader(self):
        directory = '/'.join(self.traj_loc.split('/')[0:-1])
        beef_array = pickle.load(open(directory + '/ensemble.pkl', 'r'))
        if abs(np.mean(beef_array)) < 100:
            beef_array += read(self.traj_loc).get_potential_energy()
        return beef_array       
    
    def vibs_reader(self,vib_loc):
        f = open(vib_loc+'/myjob.out','r')
        vibenergies = []
        for line in f:
            try:
                temp = line.replace('.','')
                temp = temp.replace(' ','')
                temp = temp.replace('i','')
                float(temp)
            except ValueError:
                continue
            num, meV, inv_cm, = line.split()
            if 'i' in meV:
                meV = 7
            vibenergies.append(float(meV))
        #convert from meV to eV for each mode
        vibenergies[:] = [round(ve/1000.,4) for ve in vibenergies]    
        return vibenergies

    def get_Gcorr(self,T,P=101325,verbose=False):
        if self.species_type == 'gas':
            Gcorr = self.gibbs.get_gibbs_energy(T,P,verbose=verbose)
        elif self.species_type == 'adsorbate':
            Gcorr = self.gibbs.get_helmholtz_energy(T,verbose=verbose)
        elif self.species_type == 'slab':
            Gcorr = 0
        else:
            #Think about how to throw an error here
            print "Error: ambiguous species type for function get_dGcorr.  Should be 'gas' or 'adsorbate'."
            exit()
        return Gcorr

    def get_Hcorr(self,T,verbose=False):
        if self.species_type == 'gas':
            return self.gibbs.get_enthalpy(T, verbose=verbose)
        elif self.species_type == 'adsorbate':
            return self.gibbs.get_internal_energy(T, verbose=verbose)


class AtomicConfigs:
    """
    """

    def __init__(self, dftclasses):
        self.data = dftclasses

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

        return dftclasses(out)

    def get_property(self, ppt):
        """
        Return list of property for each chargexfer
        """
        out = []
        for c in self.data:
            eval('out.append(c.%s)' % ppt)
        return np.array(out)
