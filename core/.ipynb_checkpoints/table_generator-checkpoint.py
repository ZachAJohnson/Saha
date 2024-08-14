# Zach Johnson 6/27/24
# Generates tables of Saha ionization

import numpy as np
from saha.core.saha import calculate_ionization_fractions
from hnc.hnc.constants import *

import numpy as np
from time import time

class saha_table:

	def __init__(self, plasma, n0_AU_range, T_AU_range, Nn=100, NT=100, N_ions=None):
	    self.plasma = plasma
	    self.n0_range = tuple(n0_AU_range)
	    self.T_range = tuple(T_AU_range)
	    self.Nn = Nn
	    self.NT = NT
	    self.N_ions = N_ions

	    self.gen_nT_mesh()
	    self.solve_saha()
	    self.save_mesh()

	def gen_nT_mesh(self):
	    self.n0_array = np.geomspace(*self.n0_range, num=self.Nn)
	    self.T_array = np.geomspace(*self.T_range, num=self.NT)

	def solve_saha(self):
		# Initialize empty arrays to store results
		self.Zbar_list, self.xi_list, self.χ_list = [], [], []
		self.n0_list, self.T_list = [], [] # brain hurty, easy peasy

		t0 = time()
		iters = 0
		# Iterate over the mesh grid
		for i, n0 in enumerate(self.n0_array):
			for j, T in enumerate(self.T_array):
				Zbar, xi_array, χ_array = calculate_ionization_fractions(self.plasma, n0, T, IPD=True, N_ions=self.N_ions)
				self.Zbar_list.append(Zbar)
				self.xi_list.append(xi_array)
				self.χ_list.append(χ_array)
				self.n0_list.append(n0)
				self.T_list.append(T)
				print(f"{iters/self.Nn/self.NT*100:0.3f}%", end='\r')
				iters+=1
		print(f"\n Time to finish mesh: {time()-t0:0.3e} [s], or {(time()-t0)/self.Nn/self.NT:0.3e} [s] per point ")

	def save_mesh(self):

	    mesh_save_data = np.array([np.array(self.n0_list)*AU_to_invcc, np.array(self.T_list)*AU_to_K, self.Zbar_list, *(np.array(self.xi_list).T), *(np.array(self.χ_list).T*AU_to_eV)]).T
	 
	    col_names = f'   {"n[1/cc]":12} {"T[K]":12} {"Zbar":12} '
	    for i in range(len(self.xi_list[0])):
	    	xi_name = f"x_{i}[eV]"
	    	col_names += f"{xi_name:12} "
	    for i in range(len(self.xi_list[0])):
	    	χi_name = f"χ_{i}[eV]"
	    	col_names += f"{χi_name:12} "
	    
	    header = ("# Saha solution with Stewart-Pyat α=1 IPD\n" + col_names)
	    savename = f"{self.plasma.name}_Saha.txt"
	    np.savetxt("../data/" + savename, mesh_save_data, header=header, fmt="%12.3e", comments='')

