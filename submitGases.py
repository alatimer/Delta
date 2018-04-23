#!/usr/bin/env python

import os
import glob

#gases = ['CH4','O2','H2O']
psps = ['esp','gbrv']
xcs = ['RPBE','BEEF','PBE']
pws = ['400','450','500','550','600','650','700']
#test
#pws = ['400']
#xcs = ['RPBE']
gases = ['N2O']


for gas in gases:
    for psp in psps:
        for xc in xcs:
            for pw in pws:
               traj_loc='/home/alatimer/work_dir/gases/%s/%s/%s/%s'%(gas,psp,xc,pw)
               print traj_loc
               os.system('mkdir -p %s'%(traj_loc))
               os.system('cp gas_opt.py %s'%(traj_loc))
               os.system('cp %s.traj %s/init.traj'%(gas,traj_loc))
               os.chdir(traj_loc)
               #os.system(' rm -r vib* myjob* calcdir* qn.log qn.traj uniq* node*')

               print os.getcwd()
               #for pkl in glob.glob('vib.???.pckl'):
                #   if os.stat(pkl).st_size == 0:
                 #      os.system('echo rm %s'%(pkl))
               os.system(' sbatch --job-name=$(pwd) gas_opt.py %s %s /home/alatimer/bin/psp/%s/'%(xc, pw, psp))
               os.chdir('/home/alatimer/src/Delta/')
