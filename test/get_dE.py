#!/usr/bin/env python

import pickle
from Delta import Reactant,Reaction
import sys

#H2,H2O ref
refs = {'H':{'H2':0.5},'O':{'H2O':1,'H2':-1},'C':{'CH4':1,'H2':-2}}

#H2,O2 ref
#refs = {'H':{'H2':0.5},'O':{'O2':0.5},'C':{'CH4':1,'H2':-2}}

IS = Reactant.Reactant(
    traj_loc=sys.argv[1],
    species_type = 'adsorbate',
              )

FS = Reactant.Reactant(
    traj_loc=sys.argv[2],
    species_type = 'adsorbate',
              )

rxn = Reaction.Reaction(FS=FS,IS=IS,refs=refs)
print refs
print rxn.get_dE(verbose=True)

