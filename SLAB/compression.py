# Generate a system with randomly placed polymers.
import hoomd
import hoomd.md
from hoomd import deprecated
import math

hoomd.context.initialize()
system = hoomd.init.read_gsd('init.gsd', frame=-1)
initial_L = 60;
final_L   = 30;
shrink_L = hoomd.variant.linear_interp([(0, initial_L), (5000, final_L)]) 
resize = hoomd.update.box_resize(L=shrink_L, period=1, scale_particles=True);
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=3.0, nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

all = hoomd.group.all()
hoomd.md.integrate.mode_standard(dt=0.005)
kT=4.0
integrator=hoomd.md.integrate.langevin(group=all, kT=kT, seed=42)
integrator.set_gamma('A', gamma=0.1)

outf  = "log-output_compress.log"
trajf = "compress.gsd"

hoomd.analyze.log(filename=outf, quantities=
['potential_energy',  'pressure', 'pressure_xx', 'pressure_yy', 'pressure_zz'], period=10, overwrite=True)
hoomd.dump.gsd(filename=trajf, period=1000, group=hoomd.group.all(), overwrite=True)

hoomd.run(20000)
