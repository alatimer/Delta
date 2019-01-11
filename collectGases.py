#!/usr/bin/env python

import pickle
from Delta import Reactant
from ase.io import read,write
import os

gases = ['N2O','CH4','O2','H2O','H2','CO2','CO','CH2O','CH3OH','N2','NH3','O3']
psps = ['esp','gbrv']
xcs = ['RPBE','BEEF','PBE']
#pws = ['400','550']
pws = ['400','450','500','550','600','650','700']
gc_list = []

gas_params = {'H2':[0,2,'linear'],
            'O2':[2,2,'linear'],
            'H2O':[0,2,'nonlinear'],
            'CH4':[0,12,'nonlinear'],
            'CO2':[0,2,'linear'],
            'CO':[2,2,'linear'],
            'CH2O':[0,2,'nonlinear'],
            'CH3OH':[0,1,'nonlinear'],
            'N2':[0,2,'linear'],
            'N2O':[0,2,'linear'],
            'NH3':[0,3,'nonlinear'],
            'O3':[0,2,'nonlinear'],
            }

home = os.getenv('HOME')

for psp in psps:
    for xc in xcs:
        for pw in pws:
            for gas in gases:
                directory=home+'/work_dir/gases/%s/%s/%s/%s/'%(gas,psp,xc,pw)
                traj_loc=directory+'/qn.traj'
                #if gas == 'O3':
                #    vib_loc = home+'/work_dir/gases/%s/%s/%s/%s/'%(gas,psp,xc,pw)
                #else:
                #    vib_loc=directory
                vib_loc = directory
                if os.path.exists(traj_loc)==True  and os.stat(traj_loc).st_size > 0:
                    energy = read(traj_loc).get_potential_energy()
                    gc = Reactant.Reactant(
                            traj_loc=traj_loc,
                            vib_loc=vib_loc,
                           # beef_loc=traj_loc,
                            species_type='gas',
                            surf_name=gas,
                            symmetrynumber=gas_params[gas][1],
                            geometry=gas_params[gas][2],
                            spin=gas_params[gas][0],
                            )
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
#for gc in gco.data:
#    print gc.calc_params,gc.traj_loc
