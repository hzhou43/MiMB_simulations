import hoomd
import gsd.hoomd
from hoomd import md
from hoomd import deprecated
from gsd import fl
from gsd import hoomd as gsdhoomd
from freud.box import Box
import numpy as np
import math

hoomd.context.initialize()
infile = "compress.gsd"
traj= gsdhoomd.HOOMDTrajectory(fl.GSDFile(infile, 'rb'))
box = traj[-1].configuration.box
box2 = traj[-1].configuration.box
fbox = Box(Lx=box[0], Ly=box[1], Lz=box[2])
positions_unwrap=[]
t=gsd.hoomd.open(name='elong.gsd', mode='wb')
npart=traj[-1].particles.N

count=0;
for f in traj:
	count=count+1
		
	if count==19:
		f.configuration.box[2]=150
		t.append(f)
		break
	
	if count==18:
		for j in range(f.particles.N):
			pos_unwrap = np.copy(f.particles.position)
			image = np.copy(f.particles.image)
			fbox.unwrap(pos_unwrap, image)
			positions_unwrap.append(pos_unwrap)
					
			if ((f.particles.position[j][2])>=15):
				
					f.particles.position[j][2] = f.particles.position[j][2] - 30
	
		t.append(f)

t.append(f)
