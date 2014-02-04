from mdconfig import *
import matplotlib.pylab as pylab
from math import pi

class MD_statistics:
	def __init__(self, md):
		self.md = md

	def get_total_unit_cells(self):
		total_unit_cells_x = self.md.unit_cells_x*self.md.nodes_x
		total_unit_cells_y = self.md.unit_cells_y*self.md.nodes_y
		total_unit_cells_z = self.md.unit_cells_z*self.md.nodes_z

		return [total_unit_cells_x, total_unit_cells_y, total_unit_cells_z]

	def calculate_system_length(self):
		total_unit_cells = self.get_total_unit_cells()
		system_length_x = total_unit_cells[0]*self.md.FCC_b
		system_length_y = total_unit_cells[1]*self.md.FCC_b
		system_length_z = total_unit_cells[2]*self.md.FCC_b
		return [system_length_x, system_length_y, system_length_z]

	def get_volume(self, path = "./"):
		volume = pylab.loadtxt(path+'/volume.txt')
		return volume

	def get_num_free_atoms(self, path = "./"):
		num_free_atoms = pylab.loadtxt(path+'/number_of_free_atoms.txt')
		return num_free_atoms

	def get_density(self, path = "./"):
		volume = self.get_volume(path=path)
		num_free_atoms = self.get_num_free_atoms(path=path)
		density = num_free_atoms / volume
		return density

	def get_number_flow_rate(self, path = "./"):
		filename = path+"/statistics/count_periodic.txt"
		count_periodic = pylab.loadtxt(filename)
		last_count_periodic = count_periodic[-1]
		time = last_count_periodic[0]
		count_z = last_count_periodic[3]
		number_flow_rate = count_z / time
		return number_flow_rate

	def update_new_volume(self, volume, path = "./"):
		volume_file = open(path+'/volume.txt', 'w')
		volume_file.write(str(volume))
		volume_file.close()

	def get_ideal_gas_pressure(self, temperature, path="./"):
		density = self.get_density(path=path)
		pressure = density*temperature
		return pressure

	def calculate_permeability(self, porosity, path = "./"):
		volume = self.get_volume(path=path)
		num_free_atoms = self.get_num_free_atoms(path=path)
		number_flow_rate = self.get_number_flow_rate(path=path)
		density = num_free_atoms / volume
		volume_per_atom = volume / num_free_atoms
		volumetric_flow_rate = number_flow_rate*volume_per_atom
		sigma = 3.405
		viscosity = 5.0/(16.0*sigma**2)*sqrt(self.md.mass*self.md.temperature/pi)
		
