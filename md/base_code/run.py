# This example file shows how to melt an SiO2 system

from mdconfig import *
from md_geometry import *
# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = MD()
md = program.compile(skip_compile=False, name="md")
geometry = MD_geometry(program)

if True:
	program.reset()
	program.nodes_x = 1
	program.nodes_y = 1
	program.nodes_z = 1
	program.unit_cells_x = 50
	program.unit_cells_y = 50
	program.unit_cells_z = 50

	program.prepare_new_system()
	program.run(save_state_path="states/00_initial_state")

	#program.prepare_thermostat(temperature=300, timesteps=2000, run=True, save_state_path="states/01_T_300K")
	#program.prepare_thermalize(timesteps=1000, run=True, save_state_path="states/02_thermalized")

	geometry.create_cylinders(radius=0.01, num_cylinders_per_dimension=4)

if False:
	program.load_state(path="states/03_cylinder")
	program.gravity_force = 0.001
	program.gravity_direction = 2
	program.prepare_thermalize(timesteps=5000, run=True, save_state_path="states/04_cylinder_thermalized")

	program.prepare_thermalize(timesteps=10000, run=True, save_state_path="states/05_cylinder_thermalized_more")

if True:
	program.create_movie_files = True
	program.prepare_thermalize(timesteps=100, run=True)
	program.create_movie(frames=100)