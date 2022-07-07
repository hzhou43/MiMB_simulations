# Generate a system with randomly placed polymers.
import hoomd
import hoomd.md
import gsd.hoomd
from gsd import fl
from gsd import hoomd as gsdhoomd
from hoomd import deprecated
import math
import numpy as np

hoomd.context.initialize()
system = hoomd.init.read_gsd('elong.gsd', frame=-1)

nl=hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=3.0, nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.set_params(mode="shift")

all = hoomd.group.all()
hoomd.md.integrate.mode_standard(dt=0.005)
kT = float(hoomd.option.get_user()[0])
integrator=hoomd.md.integrate.langevin(group=all, kT=kT, seed=42)
integrator.set_gamma('A', gamma=0.1)

outf  =  "log-output_equil_T" + str(kT) + ".log"
trajf =  "equil_T" + str(kT) + ".gsd"

hoomd.analyze.log(filename=outf, quantities=['potential_energy',  'pressure', 'pressure_xx', 'pressure_yy', 'pressure_zz'], period=10, overwrite=True)
gsd = hoomd.dump.gsd(filename=trajf, group=hoomd.group.all(), period=1000, overwrite=True)

hoomd.run(1e6)
