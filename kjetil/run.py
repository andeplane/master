from something import *
import numpy as np

V2 = VectorFunctionSpace(mesh, "CG", 1)
Q2 = FunctionSpace(mesh, "CG", 1)

filename_velocity = "output/v.pvd"
filename_pressure = "output/p.pvd"
file_v = File(filename_velocity)
file_p = File(filename_pressure)


def run():
	global file_v
	global file_p
	
	t,u,p = step()
	p_array = interpolate(p,Q2)

	start_time = time.time()
	file_v << u
	file_p << p_array
	elapsed_time = time.time() - start_time
	print "Files saved in ",elapsed_time, " seconds"

	print "t=",t

global boundaries
plot(boundaries)
interactive()

#while True:
#	run()