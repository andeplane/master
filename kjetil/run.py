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
	
	#u_array = interpolate(u,V2)
	p_array = interpolate(p,Q2)

	file_v << u
	file_p << p_array

	print "t=",t
	#print "v_max=",max(u_array)

while True:
	run()

plt.show()