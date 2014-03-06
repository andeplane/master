from mdconfig import *
from md_geometry import *
from md_unit_converter import *
from md_statistics import *

program = MD(dt=0.02, name="main")
md = program.compile(skip_compile=True)
geometry = MD_geometry(program)
uc = MD_unit_converter(program)
md_statistics = MD_statistics(program)
program.nodes_x = 2
program.nodes_y = 2
program.nodes_z = 2
program.unit_cells_x = 10
program.unit_cells_y = 10
program.unit_cells_z = 10
program.max_number_of_atoms = 1000000
program.max_number_of_cells = 10000
skip_start = False
if not skip_start:
	program.reset()

	program.prepare_new_system()
	program.run()

	#program.prepare_thermostat(temperature=300, timesteps=2000, run=True, save_state_path="states/01_T_300K")
	program.load_state(path="states/01_T_300K")
	program.prepare_thermalize(timesteps=10000, run=True, save_state_path="states/02_thermalized")

if not skip_start:
	program.load_state(path="states/02_thermalized")
	geometry.create_cylinders(radius=0.45, num_cylinders_per_dimension=1)
	program.save_state(path="states/03_cylinder")

if not skip_start:
	program.load_state(path="states/03_cylinder")
	program.reduce_density(relative_density = 0.01)
	program.save_state(path="states/04_cylinder_reduced_density")

print "Density: ", uc.number_density_to_si(md_statistics.get_density())

program.load_state(path="states/04_cylinder_reduced_density")
#program.prepare_thermostat(temperature=300, timesteps=2000, run=True)
program.prepare_thermalize(timesteps=1000000, run=True)

exit()
if False:
	program.load_state(path="states/04_cylinder_reduced_density")
	ideal_gas_pressure = md_statistics.get_ideal_gas_pressure(temperature=300)
	pressure_difference = 0.1*ideal_gas_pressure
	system_size = md_statistics.calculate_system_length()

	gravity_force = uc.pressure_difference_to_gravity(delta_p=pressure_difference, length=system_size[2])
	print "Gravity force, SI:", uc.acceleration_to_si(gravity_force)
	print "Gravity force, MD:", gravity_force
	program.gravity_force = gravity_force
	program.gravity_direction = 2
	program.thermostat_relaxation_time = 0.1
	program.temperature = 300
	program.thermostat_frozen_enabled = True
	program.create_movie_files = True
	program.prepare_thermalize(timesteps=10000, run=True, save_state_path="states/05_cylinder_thermalized")
	program.create_movie(frames=10000)

if True:
	#program.load_state(path="states/04_cylinder_thermalized")
	program.load_state(path="states/04_cylinder_reduced_density")
	program.create_movie_files = True
	program.prepare_thermalize(timesteps=10000, run=True)
	program.create_movie(frames=10000)