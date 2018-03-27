#!/usr/bin/env python

import pickle
from Delta import Gas,Reactant
from ase.io import read,write
import os

gases = ['CH4','O2','H2O']
psps = ['esp','gbrv']
xcs = ['RPBE','BEEF','PBE']
#pws = ['400','550']
pws = ['400','450','500','550','600','650','700']
gc_list = []

for psp in psps:
    for xc in xcs:
        for pw in pws:
            for gas in gases:
                directory='/home/alatimer/work_dir/gases/%s/%s/%s/%s/'%(gas,psp,xc,pw)
                traj_loc=directory+'/pre.traj'
                if os.path.exists(traj_loc)==True  and os.stat(traj_loc).st_size > 0:
                    energy = read(traj_loc).get_potential_energy()
                    print traj_loc
                    gc = Reactant.Reactant(
                            traj_loc=traj_loc,
                          #  vib_loc=directory,
                            species_type='gas',
                            surf_name=gas)
                    calc_check = gc.get_calc_params()
                    ##Always giving same psp tag
                    if calc_check['pw'] == pw and calc_check['xc'] == xc:
                        #way to check psps?
                        #psp_dict = calc_check['psp']
                        gc_list.append(gc)
                    else:
                        print "Warning: Calculation parameters do not match for %s with params: %s, %s, %s. Skipping."%(gas,psp,xc,str(pw))
                        print calc_check
                else:
                   traj_loc='/home/alatimer/work_dir/gases/%s/%s/%s/%s/'%(gas,psp,xc,pw)
                   print "Warning: No calculation found for %s with params: %s, %s, %s"%(gas,psp,xc,str(pw))

#remove duplicates?
#gc_list = set(gc_list)
gcs = Reactant.Reactants(gc_list)

### make pickle file
pickle.dump( gcs, open( "gases.pkl", "wb" ) )

#check that pickle file loads correctly
gco = pickle.load(open('gases.pkl','rb'))
for gc in gco.data:
    print gc.calc_params,gc.traj_loc
