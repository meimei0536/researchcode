#!/usr/bin/env python                                                                                   

from ase import io
from ase.lattice import bulk
from ase.lattice.surface import surface
from ase.io.trajectory import Trajectory
from ase.lattice.spacegroup import crystal
from ase.optimize import QuasiNewton
import numpy as np
from ase.io.trajectory import PickleTrajectory
from ase.build import cut
#import cut

a = 8.171

Bulk = crystal(('Mg','Al','O'),
               basis=[(0,0,0),
                      (0.625,0.625,0.625),
                      (0.375,0.375,0.375)],
            spacegroup = 227,
            cellpar = [a, a, a, 90, 90, 90])
Bulk.set_initial_magnetic_moments([0.0 for atom in Bulk])

#atoms = Bulk.repeat((2,2,2))
#atoms.rattle(stdev=0.001, seed=42)
io.write('1_1_1.traj',Bulk)

#al111 = surface(Bulk, (1,1,1),4, vacuum = 7.5)

#al111 = cut(Bulk,(1,0,0),(0,1,0),(0,0,1),nlayers=3)

#io.write('surface_111.traj',al111)
