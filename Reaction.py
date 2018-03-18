#Embedded file name: /scratch/users/alatimer/meoh-vs-methane/analysis/catclass.py
from ase import Atoms
from ase.io import read, write
import pickle
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from AtomicConfig import AtomicConfig
from ase.atoms import string2symbols

class Reaction:
    """
    """

    def __init__(self, 
            IS = None, #AtomicConfig type
            FS =  None, #AtomicConfig type
            refs = {'H':{'H2':0.5},'O':{'O2':0.5},'C':{'CH4':1,'H2':-2}},
            name = None, #string
            tag = '',
            ):

        self.IS = IS
        self.FS = FS
        self.name = name
        self.tag = tag
        
        if IS.calc_params == FS.calc_params:
            self.calc_params = IS.calc_params
            print self.calc_params
        else:
            print "IS and FS calc params are different."        

        self.refs = refs
        self.gases = {}
        #find somewhere else to put this
        gas_params = {'H2':[0,2,'linear'],'O2':[2,2,'linear'],'H2O':[0,2,'nonlinear'],'CH4':[0,12,'nonlinear']}
        for ref in refs:
            for gas in refs[ref]:
                if gas not in self.gases:
                    ac = AtomicConfig(
                        species_name = gas,
                        spin=gas_params[gas][0],
                        symmetrynumber = gas_params[gas][1],
                        geometry = gas_params[gas][2],
                        traj_loc='/home/alatimer/work_dir/gases/%s/%s/%s/%s/qn.traj'\
                                %(gas,self.calc_params['psp'],self.calc_params['xc'],self.calc_params['pw']),
                        beef_loc='/home/alatimer/work_dir/gases/%s/%s/%s/%s/qn.traj'\
                                %(gas,self.calc_params['psp'],self.calc_params['xc'],self.calc_params['pw']),
                        vib_loc ='/home/alatimer/work_dir/gases/%s/vibs/'%(gas),
                        species_type = 'gas',
                        calc_params = self.calc_params,
                          )
                    self.gases[gas]=ac


    def E_fun(self,ac,T,P):
        return ac.atoms.get_potential_energy()
    def G_fun(self,ac,T,P):
        return (ac.atoms.get_potential_energy()+ac.get_Gcorr(T,P))
    def get_dX(self,engfun=E_fun,T=None,P=None,verbose=False):
        comp_IS = [atom.symbol for atom in self.IS.atoms]
        comp_FS = [atom.symbol for atom in self.FS.atoms]
        dE = engfun(self.IS,T,P) - engfun(self.FS,T,P)
        
        #elements missing from FS
        if len(comp_IS) > len(comp_FS):
            comp_diff = [item for item in comp_IS if item not in comp_FS]
            for elem in comp_diff:
                for gas in self.refs[elem]:
                    dE -= self.refs[elem][gas] * engfun(self.gases[gas],T,P)
        
        #elements missing from IS
        else:
            comp_diff = [item for item in comp_FS if item not in comp_IS]
            for elem in comp_diff:
                for gas in self.refs[elem]:
                    dE += self.refs[elem][gas] * engfun(self.gases[gas],T,P)
        if verbose:
            print " Elemental Difference: ",comp_diff
        return dE
    
    def get_dE(self,verbose=False):
        return self.get_dX(engfun=self.E_fun,verbose=verbose)
    def get_dG(self,T,P=101325):
        return self.get_dX(engfun=self.G_fun,T=T,P=P)
    def get_dH(self):
        return

class Reactions:
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
