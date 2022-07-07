import hoomd
import gsd.hoomd
from gsd import hoomd as gsdhoomd
from gsd import fl
from freud.box import Box
import numpy as np
import math

hoomd.context.initialize()
kT = float(hoomd.option.get_user()[0])
infile = "equil_T" + str(kT) + ".gsd"
traj = gsdhoomd.HOOMDTrajectory(fl.GSDFile(infile, 'rb'))

num=0
count=0
outp = "Tcoor.txt"
fp=open(outp,"w+")
for f in traj:
	count=count+1
	if ((count>=500)):
		num=num+1			
		fp.write("%s\n" % count)
		for p in range(f.particles.N):
			x1 = f.particles.position[p][0] 
			x2 = f.particles.position[p][1]
			x3 = f.particles.position[p][2] 
			x4 = f.particles.typeid[p]		
			fp.write("%f %f %f %d\n" % (x1,x2,x3,x4))
	
fp.close()
