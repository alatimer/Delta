#!/usr/bin/env python

import pickle
from Delta import Gas,Reactant

#from ase import Atoms
from ase.io import read,write
import os
#import numpy as np

gases = ['CH4','H2','O2','H2O']
psps = ['esp','gbrv']
xcs = ['RPBE','BEEF','PBE']
#xcs = ['rpbe','beef','pbe']
pws = ['400','550']
#pws = ['450','500','600','650','700']
gc_list = []

for psp in psps:
    for xc in xcs:
        for pw in pws:
            for gas in gases:
                
                directory='/home/alatimer/work_dir/gases/%s/%s/%s/%s/'%(gas,psp,xc,pw)
                traj_loc=directory+'/pre.traj'
                if os.path.exists(traj_loc)==True  and os.stat(traj_loc).st_size > 0:
                #if 1==0:
                    energy = read(traj_loc).get_potential_energy()
                   # print gas,psp,xc,pw
                    print traj_loc
                    gc = Reactant.Reactant(
                            traj_loc=traj_loc,
                          #  vib_loc=directory,
                            species_type='gas',
                            surf_name=gas)
                    calc_check = gc.get_calc_params()
                    ##Always giving same psp tag
                    if calc_check['pw'] != pw or calc_check['xc'] != xc:
                        print "Warning: Calculation parameters do not match for %s with params: %s, %s, %s. Skipping."%(gas,psp,xc,str(pw))
                        print calc_check
                    else:
                        psp_dict = calc_check['psp']
                        #gc = Gas.GasCalc(name=gas,pw=calc_check['pw'],xc=calc_check['xc'],psp_dict=psp_dict,energy=energy,traj=traj,traj_loc=traj_loc)
                        
                        gc_list.append(gc)
                else:
                   traj_loc='/home/alatimer/work_dir/gases/%s/%s/%s/%s/'%(gas,psp,xc,pw)
                   print "Warning: No calculation found for %s with params: %s, %s, %s"%(gas,psp,xc,str(pw))
                   if 1==0:
                       os.system('mkdir -p %s'%(traj_loc))
                       os.system('cp gas_opt.py %s'%(traj_loc))
                       #os.system('cp %s.traj %s/init.traj'%(gas,traj_loc))
                       ##print traj_loc
                       os.chdir(traj_loc)
                       print os.getcwd()
                       os.system('sbatch --job-name=$(pwd) gas_opt.py %s %s /home/alatimer/bin/psp/%s/'%(xc, pw, psp))
                       os.chdir('/home/alatimer/src/Delta/')
                     #  print os.getcwd()

#remove duplicates
gc_list = set(gc_list)
gcs = Reactant.Reactants(gc_list)

### make pickle file
pickle.dump( gcs, open( "gases.pkl", "wb" ) )
gco = pickle.load(open('gases.pkl','rb'))

for gc in gco.data:
    print gc.calc_params,gc.traj_loc
