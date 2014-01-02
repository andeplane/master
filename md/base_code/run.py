# This example file shows how to melt an SiO2 system

from mdconfig import *

# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = MD()
md = program.compile(skip_compile=False, name="md")

program.reset()
program.nodes_x = 2
program.nodes_y = 2
program.nodes_z = 2
program.unit_cells_x = 5
program.unit_cells_y = 5
program.unit_cells_z = 5

program.prepare_new_system()
program.run(md)

program.prepare_thermostat(temperature=300, timesteps=2000, run=True)
program.prepare_thermalize(timesteps=1000, run=True)

### Create cylinder
system_length = 57.2
system_length_half = system_length/2.0
radius = system_length_half*0.9

program.create_cylinder(radius, 2, system_length_half, system_length_half, system_length_half)
### End create cylinder

program.create_movie_files = True
program.prepare_thermalize(timesteps=1000, run=True)
program.create_movie(frames=1000)
