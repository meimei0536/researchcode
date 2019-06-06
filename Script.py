#!/usr/bin/env python

import os
from ase import io, units
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength
from ase.calculators.cp2k import CP2K

inp ='''&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    CHARGE 0
    MULTIPLICITY 0
    LSD TRUE
    SURFACE_DIPOLE_CORRECTION TRUE
    SURF_DIP_DIR Z
    &SCF
       MAX_SCF  100
       EPS_SCF   1E-06
       SCF_GUESS  RESTART
       &OT  T
         MINIMIZER CG
         PRECONDITIONER  FULL_SINGLE_INVERSE
       &END OT
       &OUTER_SCF  T
         EPS_SCF    1E-06
         MAX_SCF  30
       &END OUTER_SCF
     &END SCF
     &MGRID
       NGRIDS 4
       CUTOFF  480
       REL_CUTOFF  60
     &END MGRID
     &XC
       DENSITY_CUTOFF     1.0000000000000000E-10
       GRADIENT_CUTOFF     1.0000000000000000E-10
       TAU_CUTOFF     1.0000000000000000E-10
       &XC_GRID
         XC_SMOOTH_RHO  NN50
         XC_DERIV  SPLINE3
       &END XC_GRID
       &XC_FUNCTIONAL
       &END XC_FUNCTIONAL
       &VDW_POTENTIAL
         DISPERSION_FUNCTIONAL PAIR_POTENTIAL
         &PAIR_POTENTIAL
           TYPE DFTD3
           PARAMETER_FILE_NAME dftd3.dat
           REFERENCE_FUNCTIONAL PBE
         &END PAIR_POTENTIAL
       &END VDW_POTENTIAL
     &END XC
  &END DFT
  &SUBSYS
    &KIND Mg
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE
    &END KIND
    &KIND Al
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE
    &END KIND
    &KIND O
      BASIS_SET TZV2P-MOLOPT-GTH
      POTENTIAL GTH-PBE
    &END KIND
    &KIND Ir
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE
    &END KIND
    &KIND C
      BASIS_SET TZV2P-MOLOPT-GTH
      POTENTIAL GTH-PBE
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
'''

if os.path.exists('qn.traj') and os.path.getsize('qn.traj')!=0:
    atoms=io.read('qn.traj',index=-1)
else:
    atoms = io.read('InitialGeom.xyz')

cell=([18.21500974,0,0],[-7.2860039,8.92349591,0],[0,0, 34.94966576])
atoms.set_cell(cell)
atoms.set_pbc([1,1,1])

constraint_s = " ".join([str(i+1) for i in range(len(atoms)) if atoms[i].position[2]<9])
print constraint_s
constraint = FixAtoms(indices = [atom.index for atom in atoms if atom.position[2]<9])
c = FixBondLength(274,276)

atoms.set_constraint([constraint,c])
io.write('InitialGeom.traj',atoms)

SYSNAME = 'MgAl2O4'

calc = CP2K(label = SYSNAME,
            xc = 'PBE',
            cutoff = None,
            basis_set_file = None,
            potential_file= None,
            basis_set = None,
            pseudo_potential = None,
            stress_tensor= False,
            max_scf = None,
            inp = inp,
            print_level='LOW')

atoms.set_calculator(calc)
from ase.io.trajectory import Trajectory
Traj = Trajectory('qn.traj','a',atoms)
dyn = QuasiNewton(atoms,trajectory=Traj,logfile='qn.log',restart='qn.pckl')
dyn.run(fmax=0.1)
