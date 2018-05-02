#!/usr/bin/env python

import pickle
from Delta import Reactant,Reaction
import sys
import os

#H2,H2O ref
refs = {'H':{'H2':0.5},'O':{'H2O':1,'H2':-1},'C':{'CH4':1,'H2':-2}}

#H2,O2 ref
#refs = {'H':{'H2':0.5},'O':{'O2':0.5},'C':{'CH4':1,'H2':-2}}

home = os.getenv('HOME')

IS = Reactant.Reactant(
    traj_loc=home+'/src/Delta/test/COOH/qn.traj',
    vib_loc=home+'/src/Delta/test/COOH/',
    species_type = 'adsorbate',
              )

FS = Reactant.Reactant(
    traj_loc=home+'/src/Delta/test/slab/qn.traj',
    species_type = 'slab',
              )

rxn = Reaction.Reaction(FSs=[FS],ISs=[IS],refs=refs)
print refs
print "dE: ",rxn.get_dE(verbose=False)
print "dG: ",rxn.get_dG(T=300,P=101325,verbose=False)
