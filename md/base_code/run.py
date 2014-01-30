# This example file shows how to melt an SiO2 system

from mdconfig import *
from md_geometry import *
from md_unit_converter import *
from md_statistics import *
# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = MD(dt=0.01)
md = program.compile(skip_compile=False, name="md")
geometry = MD_geometry(program)
uc = MD_unit_converter(program)
md_statistics = MD_statistics(program)
program.nodes_x = 1
program.nodes_y = 1
program.nodes_z = 1
program.unit_cells_x = 10
program.unit_cells_y = 10
program.unit_cells_z = 10

if True:
	program.reset()

	program.prepare_new_system()
	program.run()

	program.prepare_thermostat(temperature=300, timesteps=2000, run=True, save_state_path="states/01_T_300K")
	program.prepare_thermalize(timesteps=1000, run=True, save_state_path="states/02_thermalized")
	geometry.create_cylinders(radius=0.4, num_cylinders_per_dimension=1)
	program.save_state(path="states/03_cylinder")

if False:
	program.load_state(path="states/03_cylinder")
	ideal_gas_pressure = md_statistics.get_ideal_gas_pressure(temperature=300)
	pressure_difference = 0.05*ideal_gas_pressure
	system_size = md_statistics.calculate_system_length()

	gravity_force = uc.pressure_difference_to_gravity(delta_p=pressure_difference, length=system_size[2])
	print "Gravity force, SI:", uc.acceleration_to_si(gravity_force)
	print "Gravity force, MD:", gravity_force
	program.gravity_force = gravity_force
	program.gravity_direction = 2
	program.thermostat_relaxation_time = 0.1
	program.temperature = 300
	program.thermostat_frozen_enabled = True
	program.prepare_thermalize(timesteps=10000, run=True, save_state_path="states/04_cylinder_thermalized")

if True:
	#program.load_state(path="states/04_cylinder_thermalized")
	program.load_state(path="states/03_cylinder")
	program.create_movie_files = True
	program.prepare_thermalize(timesteps=1000, run=True)
	program.create_movie(frames=1000)