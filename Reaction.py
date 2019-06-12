from ase import Atoms
from ase.io import read, write
import pickle
import numpy as np
import os
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from Reactant import Reactant
from ase.atoms import string2symbols
import copy

class Reaction:
    """
    """

    def __init__(self, 
            ISs = None, #Reactant type
            FSs =  None, #Reactant type
            dft_refs = {'H':{'H2':0.5},'O':{'O2':0.5},'C':{'CH4':1.,'H2':-2.}},
            name = None, #string
            tag = '',
            electrochemical=False,
            ignore=None, ###????
            #Default calc params for gas eforms
            calc_params={
                'xc': 'BEEF', 
                'psp': {'H': 'bd57da0722753ddfd9833a65b40bb853', 
                        'C': 'ab446d787031eba755f1d5ee46644bb0', 
                        'O': 'e057bbcf3229fb5d8c36c0c92f9568b4'}, 
                'pw': '550', 
                'kpts': '111'},
            ):

        self.ISs = copy.copy(ISs)
        self.FSs = copy.copy(FSs)
        self.name = name
        self.tag = tag
        self.dft_refs = dft_refs
        self.electrochemical=electrochemical
        self.gases = {}
        self.calc_params = calc_params
        self.default_params = calc_params
        #checking calc params
        for i,state in enumerate(self.ISs+self.FSs):
            if isinstance(state,basestring):
                continue
            if i==0:
                #Set initial calc_params to be that of first state
                self.calc_params=state.calc_params
                continue
            #make sure calculation parameters are equivalent in all states
            if (    self.calc_params['pw'] == state.calc_params['pw'] and
                    self.calc_params['xc'] == state.calc_params['xc'] ):
                if (    self.calc_params['kpts'] != state.calc_params['kpts'] and
                        self.ISs[0].species_type != 'bulk'):
                    print "Warning: calc params are different. Make sure a gas is not the first state passed."    
                    print self.calc_params,'\n', state.calc_params
                    exit()

                for elem in state.calc_params['psp']:
                    if (    elem in self.calc_params['psp'] and 
                            self.calc_params['psp'][elem] != self.calc_params['psp'][elem]):
                        print 'PSPs not equivalent'
                        exit()
                    elif elem not in self.calc_params['psp']:
                        self.calc_params['psp'][elem]=state.calc_params['psp'][elem]
            else:
                    print "Warning: calc params are different. Make sure a gas is not the first state passed."    
                    print self.calc_params,'\n', state.calc_params
                    exit()
        
        #Open gas DB
        classloc =  '/'.join(__file__.split('/')[0:-1])
        gasDB = pickle.load(open(classloc+'/gases.pkl','rb'))

        #if any states are strings, convert them to relevant reactant object
        for states in [self.ISs,self.FSs]:
            for i,state in enumerate(states):
                if isinstance(state,basestring):
                    #if all states are strings, calc params are default (provided)
                    #if not all states are strings, calc params are those of non-string
                    g = gasDB.filter(lambda x: self.calc_params['xc'] == x.calc_params['xc'])
                    g = g.filter(lambda x: self.calc_params['pw'] == x.calc_params['pw'])
                    g = g.filter(lambda x: state == x.surf_name)
                    for elem in string2symbols(state):
                        if  elem in self.calc_params['psp']:
                            g = g.filter(lambda x: self.calc_params['psp'][elem] == x.calc_params['psp'][elem])
                        elif elem not in self.calc_params['psp']:
                            g = g.filter(lambda x: self.default_params['psp'][elem] == x.calc_params['psp'][elem])
                    if len(g.data)>1:
                        print "Warning: Ambiguous gas phase reference"
                    states[i]= g.data[0]


        #Find elemental difference between IS and FS
        self.comp_dict = {}
        for IS in self.ISs:
            for atom in IS.atoms:
                if self.comp_dict.has_key(atom.symbol)==True:
                    self.comp_dict[atom.symbol]+=1
                else:
                    self.comp_dict[atom.symbol]=1
        for FS in self.FSs:
            for atom in FS.atoms:
                if self.comp_dict.has_key(atom.symbol)==True:
                    self.comp_dict[atom.symbol]-=1
                else:
                    self.comp_dict[atom.symbol]=-1

        #Collecting reactants of required reference gases
        #This info will be moved to gas DB
        gas_params = {'H2':[0,2,'linear'],'O2':[2,2,'linear'],'H2O':[0,2,'nonlinear'],'CH4':[0,12,'nonlinear'],'CH3OH':[0,1,'nonlinear']}
        
        gasDB = gasDB.filter(lambda x: self.calc_params['xc'] == x.calc_params['xc'])
        gasDB = gasDB.filter(lambda x: self.calc_params['pw'] == x.calc_params['pw'])
        
        #Finding ref gases with matching psps, janky...
        i = 0
        while i<4:
          i+=1
          for elem in self.dft_refs:
            if elem in self.calc_params['psp']:
              for gas in dft_refs[elem]:
                if elem in gas:
                  g = gasDB.filter(lambda x: gas == x.surf_name)
                  g = g.filter(lambda x: self.calc_params['psp'][elem] == x.calc_params['psp'][elem])
                  if len(g.data)>1:
                      print "Warning: Unclear choice of reference psp:"
                      g.print_all()
                  self.gases[gas] = g.data[0]
                  for elem2 in string2symbols(gas):
                    if elem2 not in self.calc_params:
                      self.calc_params['psp'][elem2] = self.gases[gas].calc_params['psp'][elem2]
        return

    def E_fun(self,ac,T,P):
        #T,P ignored
        if ac !=None:
            eng = ac.atoms.get_potential_energy()
        else:
            eng = 0
        return eng
    def G_fun(self,ac,T,P):
        if ac !=None:
            eng = (ac.atoms.get_potential_energy()+ac.get_Gcorr(T,P))
        else:
            eng = 0
        return eng
    
    def H_fun(self,ac,T,P):
        #P ignored
        if ac !=None:
            eng = (ac.atoms.get_potential_energy()+ac.get_Hcorr(T))
        else:
            eng = 0
        return eng
    
 

    def get_dX(self,engfun=E_fun,T=None,P=None,verbose=False):
        dE = 0
        for i,IS in enumerate(self.ISs):
            if P!=None: #dG or dH
                dE-=engfun(IS,T,P[i])
            else: #dE
                dE-=engfun(IS,T,P)
        for i,FS in enumerate(self.FSs):
            if P!=None: #dG or dH
                dE+=engfun(FS,T,P[len(self.ISs)+i])
            else: #dE
                dE+=engfun(FS,T,P)
        for elem in self.comp_dict:
            if self.comp_dict[elem]!=0:
                #All references at P=101325 except liquid H2O if echem==True. 
                #To input other pressures, specify manually
                for gas in self.dft_refs[elem]:
                    ##for liquid water
                    if gas == 'H2O' and self.electrochemical==True:
                        dE += self.comp_dict[elem]*(self.dft_refs[elem][gas] * engfun(self.gases[gas],T,101325*0.035))
                    else:
                        dE += self.comp_dict[elem]*(self.dft_refs[elem][gas] * engfun(self.gases[gas],T,101325))
        if verbose:
            print " Elemental Difference: ",self.comp_dict
        return dE
    
    def get_dE(self,verbose=False):
        return self.get_dX(engfun=self.E_fun,verbose=verbose)
    def get_dG(self,T,P=101325,verbose=False):
        if type(P)!=list:
            P = [P]*(len(self.ISs)+len(self.FSs))
        elif len(P)!=(len(self.ISs)+len(self.FSs)):
            print "Error: not enough pressures provided for all reactants"
            exit()
        return self.get_dX(engfun=self.G_fun,T=T,P=P,verbose=verbose)
    def get_dH(self,T,P=101325):
        if type(P)!=list:
            P = [P]*(len(self.ISs)+len(self.FSs))
        elif len(P)!=(len(self.ISs)+len(self.FSs)):
            print "Error: not enough pressures provided for all reactants"
            exit()
        return self.get_dX(engfun=self.H_fun,T=T,P=P,verbose=False)

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
