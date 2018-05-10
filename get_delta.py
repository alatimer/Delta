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

COOH = Reactant.Reactant(
    species_name='COOH',
    traj_loc=home+'/src/Delta/test/COOH/qn.traj',
    vib_loc=home+'/src/Delta/test/COOH/',
    species_type = 'adsorbate',
              )

slab = Reactant.Reactant(
    species_name='slab',
    traj_loc=home+'/src/Delta/test/slab/qn.traj',
    species_type = 'slab',
              )

rxn = Reaction.Reaction(FSs=[slab],ISs=[COOH],refs=refs)
print refs
print "dE: ",rxn.get_dE(verbose=False)
print "dG: ",rxn.get_dG(T=300,P=101325,verbose=False)

print "Testing gas, explicit"

FS = Reactant.Reactant(
    species_name='CH3OH',
    traj_loc=home+'/src/Delta/test/CH3OHg/qn.traj',
    vib_loc=home+'/src/Delta/test/CH3OHg/',
    species_type = 'gas',
    spin=0,
    geometry='nonlinear',
    symmetrynumber=1,
              )

rxn = Reaction.Reaction(FSs=[FS],ISs=[],refs=refs)
print refs
print "dE: ",rxn.get_dE(verbose=False)
print "dG: ",rxn.get_dG(T=300,P=101325,verbose=False)

print "Testing gas, implicit"

rxn = Reaction.Reaction(FSs=['H2'],ISs=[],refs=refs)
print refs
print "dE: ",rxn.get_dE(verbose=False)
print "dG: ",rxn.get_dG(T=300,P=101325,verbose=False)


print "Testing diff pw"

IS = Reactant.Reactant(
    species_name='COOH',
    traj_loc=home+'/src/Delta/test/COOH-400/qn.traj',
    vib_loc=home+'/src/Delta/test/COOH-400/',
    species_type = 'adsorbate',
              )

FS = Reactant.Reactant(
    species_name='slab',
    traj_loc=home+'/src/Delta/test/slab/qn.traj',
    species_type = 'slab',
              )

#rxn = Reaction.Reaction(FSs=[FS],ISs=[IS],refs=refs)
#print refs
#print "dE: ",rxn.get_dE(verbose=False)
#print "dG: ",rxn.get_dG(T=300,P=101325,verbose=False)


print "Testing diff psps"

IS = Reactant.Reactant(
    species_name='COOH',
    traj_loc=home+'/src/Delta/test/COOH-esp/qn.traj',
    vib_loc=home+'/src/Delta/test/COOH-esp/',
    species_type = 'adsorbate',
              )

FS = Reactant.Reactant(
    species_name='slab',
    traj_loc=home+'/src/Delta/test/slab/qn.traj',
    species_type = 'slab',
              )

#rxn = Reaction.Reaction(FSs=[FS],ISs=[IS],refs=refs)
#print "dE: ",rxn.get_dE(verbose=False)
#print "dG: ",rxn.get_dG(T=300,P=101325,verbose=False)

