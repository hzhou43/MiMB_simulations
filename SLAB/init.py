import hoomd
from hoomd import md
from hoomd import deprecated
hoomd.context.initialize()

polymer1 = dict(bond_len=1.2, type=['A']*1, bond="linear", count=20000)
system = hoomd.deprecated.init.create_random_polymers(box=hoomd.data.boxdim(L=60), polymers=[polymer1], separation=dict(A=0.25), seed=52);

nl = md.nlist.cell(check_period=1);
lj = md.pair.lj(r_cut=3.0, nlist=nl);

lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0);
all = hoomd.group.all();

hoomd.md.integrate.mode_minimize_fire(group=all, dt=0.001, Nmin=5, finc=1.1, fdec=0.5,
				      alpha_start=0.1, falpha=0.99, ftol=0.1, 
				      Etol=1e-05, min_steps=10000)
hoomd.analyze.log(filename="log-output_init.log", quantities=['temperature', 'potential_energy'], period=10, overwrite=True)
hoomd.dump.gsd("init.gsd", period=1000, group=all, overwrite=True);

hoomd.run(5000)
