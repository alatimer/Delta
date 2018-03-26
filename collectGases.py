#!/usr/bin/env python

import pickle
from Delta import Gas,AtomicConfig

#from ase import Atoms
from ase.io import read,write
import os
#import numpy as np

gases = ['CH4','H2','O2','H2O']
psps = ['esp','gbrv']
xcs = ['rpbe','beef','pbe']
pws = ['400','550']

gc_list = []

for psp in psps:
    for xc in xcs:
        for pw in pws:
            for gas in gases:
                traj='/home/alatimer/work_dir/gases/%s/%s/%s/%s/qn.traj'%(gas,psp,xc,pw)
                if os.path.exists(traj)==True  and os.stat(traj).st_size > 0:
                    energy = read(traj).get_potential_energy()
                   # print gas,psp,xc,pw
                    calc_check = AtomicConfig.AtomicConfig(traj_loc=traj).get_calc_params()
                    #print calc_check
                    if calc_check['pw'] != pw or calc_check['xc'] != xc:
                        print "Warning: Calculation parameters do not match for%s with params: %s, %s, %s. Skipping."%(gas,psp,xc,str(pw))
                    else:
                        psp_dict = calc_check['psp']
                        gc = Gas.GasCalc(name=gas,pw=pw,xc=xc,psp_dict=psp_dict,energy=energy)
                        gc_list.append(gc)
                else:
                   print "Warning: No calculation found for %s with params: %s, %s, %s"%(gas,psp,xc,str(pw))

gcs = Gas.GasCalcs(gc_list)

### make pickle file
pickle.dump( gcs, open( "gases.pkl", "wb" ) )
gco = pickle.load(open('gases.pkl','rb'))

for gc in gco.data:
    print gc.name,gc.energy,gc.psp_dict
