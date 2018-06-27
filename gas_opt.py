#!/usr/bin/env python

#SBATCH -p iric
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --time=15:00:00
#SBATCH --mail-type=FAIL,END                            #get emailed about job BEGIN, END, or FAIL
#SBATCH --mail-user=allegralatimer@gmail.com
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks-per-node=16

from ase import *
from ase.vibrations import Vibrations
from ase.optimize import QuasiNewton
from espresso import espresso
from espresso.vibespresso import vibespresso
from ase.dft.bee import BEEF_Ensemble
import cPickle as pickle
from ase.io import read,write,Trajectory
import sys
import os
import pdos,bader

if len(sys.argv) < 3:
    print "USAGE: gas_opt.py xc pw psppath"
    exit()

if os.path.exists('qn.traj')==True  and os.stat("qn.traj").st_size > 0:
        atoms =read('qn.traj')
        atoms.write('pre.traj')
else:
        atoms = read('init.traj')

xc = sys.argv[1]
pw = int(sys.argv[2])
dw = pw*10
psppath = sys.argv[3]
dipole = {'status':False}
spinpol = True
output = {'removesave':True}
outdir = 'calcdir'

save_pdos_pkl = False
save_cube = False
save_cd = False

for atom in atoms:
    if atom.symbol == 'H':
        atom.magmom = 0
    else:
        atom.magmom = 1.5

calc = espresso(pw=pw,	#plane-wave cutoff
                dw=dw,		#density cutoff
                xc=xc,		#exchange-correlation functional
                kpts=(1,1,1), #k-point sampling
                nbands=-30,	#10 extra bands besides the bands needed to hold
                					#the valence electrons
                sigma=0.1,
                psppath = psppath,
                spinpol=spinpol,
                convergence= {'energy':1e-5,
					               'mixing':0.1,
					               'nmix':10,
					               'mix':4,
					               'maxsteps':500,
					               'diag':'david'
					                },	#convergence parameters
                output=output,
                outdir=outdir)	#output directory for Quantum Espresso files

if xc=='BEEF':
    calc.beefensemble = True
    calc.printensemble = True

atoms.set_calculator(calc)

dyn = QuasiNewton(atoms, logfile='qn.log', trajectory='qn.traj')
dyn.run(fmax=0.05)
calc.stop()

if xc == 'BEEF':
    ens = BEEF_Ensemble(calc)
    ens_e = ens.get_ensemble_energies()
    ens.write('ensemble.bee')
    pickle.dump(ens_e,open('ensemble.pkl','w'))

####################### PDOS ###################

pdos_traj = Trajectory('pdos_qn.traj','w',atoms)
pdos.pdos(atoms,outdir=outdir,save_pkl=save_pdos_pkl,spinpol=spinpol)
pdos_traj.write() #final image will have charges and duplicate geometry

bader_traj = Trajectory('bader_qn.traj','w',atoms)
bader.bader(atoms,outdir=outdir,save_cube=save_cube,save_cd=save_cd,spinpol=spinpol)
bader_traj.write()

################## VIBRATIONS ######################
vibrateatoms=[]
for i in range(len(atoms)):
        vibrateatoms.append(i)

# Calculate vibrations
calc2 = vibespresso(pw=pw,	#plane-wave cutoff
                dw=dw,		#density cutoff
                xc=xc,		#exchange-correlation functional
                kpts=(1,1,1), #k-point sampling
                nbands=-30,	#10 extra bands besides the bands needed to hold
                					#the valence electrons
                sigma=0.1,
                psppath = psppath,
                spinpol=spinpol,
                convergence= {'energy':1e-5,
					               'mixing':0.1,
					               'nmix':10,
					               'mix':4,
					               'maxsteps':500,
					               'diag':'david'
					                },	#convergence parameters
                output=output,
                outdirprefix='calcdirv')	#output directory for Quantum Espresso files

atoms.set_calculator(calc2)

vib = Vibrations(atoms,indices=vibrateatoms,delta=0.03)
vib.run()
vib.summary(method='standard')

# Make trajectory files to visualize the modes.
for mode in range(len(vibrateatoms)*3):
    vib.write_mode(mode)
