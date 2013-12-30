# This example file shows how to melt an SiO2 system

from mdconfig import *

# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = MD()
md = program.compile(skip_compile=False, name="job1")

program.reset()
program.nodes_x = 1
program.nodes_y = 1
program.nodes_z = 1

program.prepare_new_system()
program.run(md)

program.timesteps = 1000
#program.create_movie_files = True
program.create_config_file()
program.run(md)
program.create_movie(frames=1000)