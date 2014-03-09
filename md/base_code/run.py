from mdconfig import *
from md_geometry import *
from md_unit_converter import *
from md_statistics import *

program = MD(dt=0.02, name="main")
md = program.compile(skip_compile=True)
geometry = MDGeometry(program)
unit_converter = MDUnitConverter(program)
md_statistics = MDStatistics(program, unit_converter)
program.nodes_x = 2
program.nodes_y = 2
program.nodes_z = 2
program.unit_cells_x = 5
program.unit_cells_y = 5
program.unit_cells_z = 5
program.max_number_of_cells = 10000
T = unit_converter.temperature_from_si(300)

if True:
	program.reset()
	program.prepare_new_system()
	program.run()
	for i in range(5):
		print "### Applying thermostat, T=300K ###"
		program.prepare_thermostat(temperature=T, timesteps=2000, run=True)
		print "### Thermalizing ... ###"
		program.prepare_thermalize(timesteps=2000, run=True)

	print "### Applying thermostat, T=300K ###"
	
	program.prepare_thermostat(temperature=T, timesteps=200, run=True, save_state_path="states/01_T_300K")
	print "### Thermalizing ... ###"
	program.prepare_thermalize(timesteps=200, run=True, save_state_path="states/02_thermalized")

program.load_state(path="states/02_thermalized")
print "### Creating cylinder ###"
geometry.create_cylinders(radius=0.45, num_cylinders_per_dimension=1)
program.save_state(path="states/03_cylinder")

program.test_mode = False
initial_density = unit_converter.number_density_to_si(md_statistics.get_density(path="states/03_cylinder"))
densities = [2e25, 4e25, 6e25, 8e25, 1e26, 3e26, 5e26, 7e26, 9e26, 1e27, 5e27, 1e28]
for i in range(len(densities)):
	print "### Loading cylinder ###"
	# Reload the cylinder again
	program.load_state(path="states/03_cylinder")
	density = densities[len(densities)-1-i]
	ratio = density / initial_density
	
	print "### Reducing density to %e ###" % (density)
	program.reduce_density(relative_density = ratio)
	program.max_number_of_atoms = 1000000
	print "### Applying thermostat, 300K ###"
	for i in range(3):
		program.prepare_thermostat(temperature=T, timesteps=2000, run=True)
		program.prepare_thermalize(timesteps=2000, run=True)
	
	print "### Thermalizing ... ###"
	program.prepare_thermalize(timesteps=10000, run=True)
	
	ideal_gas_pressure = md_statistics.get_ideal_gas_pressure(temperature=T)
	print "### Ideal gas pressure= %e ###" % (unit_converter.pressure_to_si(P=ideal_gas_pressure))
	pressure_difference = 0.1*ideal_gas_pressure
	system_size = md_statistics.calculate_system_length()
	gravity_force = unit_converter.pressure_difference_to_gravity(delta_p=pressure_difference, length=system_size[2])
	
	program.gravity_force = gravity_force
	program.gravity_direction = 2
	program.thermostat_relaxation_time = 0.1
	
	state_base_name = "states/%e/" % (density)
	
	print "### Thermalizing with gravity ... ###"
	program.prepare_frozen_thermostat(temperature=T, timesteps=20000, run=True, save_state_path=state_base_name+"04_thermalized")
	print "### Sampling statistics ###"
	program.prepare_frozen_thermostat(temperature=T, timesteps=50000, run=True, save_state_path=state_base_name+"05_sampling")
