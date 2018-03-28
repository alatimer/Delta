#Embedded file name: /scratch/users/alatimer/meoh-vs-methane/analysis/catclass.py
from ase import Atoms
from ase.io import read, write
import pickle
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from Reactant import Reactant
from ase.atoms import string2symbols

class Reaction:
    """
    """

    def __init__(self, 
            IS = None, #Reactant type
            FS =  None, #Reactant type
            refs = {'H':{'H2':0.5},'O':{'O2':0.5},'C':{'CH4':1.,'H2':-2.}},
            name = None, #string
            tag = '',
            ):

        self.IS = IS
        self.FS = FS
        self.name = name
        self.tag = tag
        self.refs = refs
        self.gases = {}
       
        #make sure calculation parameters are equivalent in IS and FS
        if (    IS.calc_params['pw'] == FS.calc_params['pw'] and
                IS.calc_params['xc'] == FS.calc_params['xc'] and
                IS.calc_params['kpts'] == FS.calc_params['kpts'] ):
            #need a way to check psps too
            if IS.calc_params['psp'] > FS.calc_params['psp']:
                self.calc_params = IS.calc_params
            else:
                self.calc_params = FS.calc_params
        else:
            print "Error: IS and FS calc params are different."    
            print "IS: ", self.IS.calc_params 
            print "FS: ", self.FS.calc_params
            exit()

        #Find elemental difference between IS and FS
        self.comp_dict = {}
        for atom in self.IS.atoms:
            if self.comp_dict.has_key(atom.symbol)==True:
                self.comp_dict[atom.symbol]+=1
            else:
                self.comp_dict[atom.symbol]=1
        for atom in self.FS.atoms:
            if self.comp_dict.has_key(atom.symbol)==True:
                self.comp_dict[atom.symbol]-=1
            else:
                self.comp_dict[atom.symbol]=-1
        #print "comp diff:",self.comp_dict

        #Collecting reactants of required reference gases
        #This info will be moved to gas DB
        gas_params = {'H2':[0,2,'linear'],'O2':[2,2,'linear'],'H2O':[0,2,'nonlinear'],'CH4':[0,12,'nonlinear']}
        
        #Reference gases
        gasDB = pickle.load(open('/home/alatimer/src/Delta/gases.pkl','rb'))
        gasDB = gasDB.filter(lambda x: self.calc_params['xc'] == x.calc_params['xc'])
        gasDB = gasDB.filter(lambda x: self.calc_params['pw'] == x.calc_params['pw'])
        #print gasDB.data

        #janky...
        i = 0
        while i<4:
          i+=1
          for elem in self.refs:
            if elem in self.calc_params['psp']:
              for gas in refs[elem]:
                if elem in gas:
                  g = gasDB.filter(lambda x: gas == x.surf_name)
                  #g.print_all()
                  g = g.filter(lambda x: self.calc_params['psp'][elem] == x.calc_params['psp'][elem])
                  #print g.data
                  if len(g.data)>1:
                      print "Warning: Unclear choice of reference psp:"
                      g.print_all()
                  #print g.data
                  self.gases[gas] = g.data[0]
                  for elem2 in string2symbols(gas):
                    if elem2 not in self.calc_params:
                      self.calc_params['psp'][elem2] = self.gases[gas].calc_params['psp'][elem2]
        #print self.calc_params

    def E_fun(self,ac,T,P):
        return ac.atoms.get_potential_energy()
    def G_fun(self,ac,T,P):
        return (ac.atoms.get_potential_energy()+ac.get_Gcorr(T,P))
    def get_dX(self,engfun=E_fun,T=None,P=None,verbose=False):
        comp_IS = [atom.symbol for atom in self.IS.atoms]
        comp_FS = [atom.symbol for atom in self.FS.atoms]
        dE = engfun(self.FS,T,P) - engfun(self.IS,T,P)
        for elem in self.comp_dict:
            if self.comp_dict[elem]!=0:
                for gas in self.refs[elem]:
                    dE += self.comp_dict[elem]*(self.refs[elem][gas] * engfun(self.gases[gas],T,P))
        if verbose:
            print " Elemental Difference: ",self.comp_dict
        return dE
    
    def get_dE(self,verbose=False):
        return self.get_dX(engfun=self.E_fun,verbose=verbose)
    def get_dG(self,T,P=101325,verbose=False):
        return self.get_dX(engfun=self.G_fun,T=T,P=P,verbose=verbose)
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
