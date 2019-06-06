#!/usr/bin/env python

from ase import Atoms, Atom, io
from ase.build import fcc111
from ase.calculators.vasp import Vasp
import os
from ase.constraints import FixAtoms
from ase.build import surface
from math import sqrt
from ase.optimize import QuasiNewton
from ase.io.trajectory import Trajectory

# Surface                                                                
atoms = io.read('../../bulk/Ni_VP4/CONTCAR')
atoms.append(Atom('P',(0.018,4.010,11.755)))
atoms.center(vacuum = 7.5, axis = 2)
print atoms
print 'Unit cell array \n', atoms.get_cell()

constraint = FixAtoms(mask=[atom.position[2]-min(atoms.positions[:,2]) < 6.2 for atom in atoms])
atoms.set_constraint(constraint)

atoms.set_initial_magnetic_moments([5 if a.symbol == 'Ni' else 0.6 for a in atoms])
print atoms.get_initial_magnetic_moments()

calc = Vasp(xc = 'PBE',
            setups = {'base': 'recommended', 'Ni': '_pv'},
            gamma = True,
            kpts = (4,4,1),
            ispin = 2,
            #prec = 'Accurate',
            encut = 400,
            ismear = -5,
            sigma = 0.1,
            nelm = 100,
            ibrion = 2,
            isif = 0,
            nsw = 0,
            idipol = 3,
            ediff = 1E-5)
if os.path.exists('qn.traj') and os.path.getsize('qn.traj')!=0:
    atoms=io.read('qn.traj',index=-1)

atoms.set_calculator(calc)

io.write('InitialGeom.traj',atoms)

Traj=Trajectory('qn.traj','a',atoms)
dyn = QuasiNewton(atoms,trajectory=Traj,logfile='qn.log',restart='qn.pckl')
dyn.run(fmax=0.05)
