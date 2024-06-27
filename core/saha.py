# Zach Johnson 6/27/24
# Saha Ionization code

from hnc.hnc.constants import *
from hnc.hnc.misc import More_TF_Zbar

import numpy as np
from scipy.optimize import minimize

from pandas import read_csv

from scipy.interpolate import interp1d

class plasma():

	def __init__(self, name, Z, χ0_AU_array, g_degeneracy_array):
		self.name = name
		self.Z = Z
		self.χ0_AU_array = χ0_AU_array
		self.g_degeneracy_array = g_degeneracy_array


def λD(T):
    return np.sqrt(2*π/(m_e*T))

def ΔU_SP(ionization_fractions, nn, Ti ): # Stewart-Pyat as written by Crowley assuming α=1, so Λ=Λ_tilde in https://arxiv.org/pdf/1309.1456
    Zj = np.arange(len(ionization_fractions))
    Zp = get_Zp(ionization_fractions)
    ne = get_ne(nn, ionization_fractions)
    rj = get_rj(Zj, ne)
    Γj_array = get_Γj(Zj, Zp, rj, Ti)
    Λ_array = get_Λ(Γj_array)
    return np.nan_to_num( - Ti/(2*Zp)*( (1 + Λ_array)**(2/3) - 1 ), nan=0)

def get_Λ(Γj): # NOT Λ_tilde
    Λ = (3*Γj)**(3/2)
    return Λ
    
def get_Zp(ionization_fractions):
    Zj = np.arange(len(ionization_fractions))
    return np.sum(ionization_fractions*Zj**2)/np.sum(ionization_fractions*Zj)

def get_rj(Zj, ne): #ion sphere radius
    return ( 3*Zj/(4*π*ne) )**(1/3)

def get_Γj(Zj, Zp, rj, Ti): # weird strong coupling parameter
    Γj = Zj*Zp/(rj*Ti)
    return Γj

def get_ne(nn, ionization_fractions):
    return nn*get_Zbar(ionization_fractions)

def get_Zbar(ionization_fractions):
    Zj = np.arange(len(ionization_fractions))
    return np.sum(ionization_fractions * Zj)


def saha_equation(ne, T, Δε, degeneracy_ratio):
    # Saha equation in atomic units
    saha_ratio = degeneracy_ratio * 2*λD(T)**-3 * np.exp(- Δε / T) / ne
    test = saha_ratio
    return test#saha_ratio

def get_ionization_fractions(plasma, ionization_fractions, nn, T, IPD = True):
    Z, χ0_AU_array, g_degeneracy_array = plasma.Z, plasma.χ0_AU_array, plasma.g_degeneracy_array

    N_atoms = len(ionization_fractions)
    N_ions = N_atoms-1

    unnormalized_fractions = [1]
    ne = get_ne(nn, ionization_fractions)
    for i in range(N_ions):
        degeneracy_ratio = g_degeneracy_array[i + 1] / g_degeneracy_array[i]
        if IPD==True:
            ipd_energies = ΔU_SP(ionization_fractions, nn, T )
        else:
            ipd_energies = np.zeros_like(ionization_fractions)
        χ_AU_array = χ0_AU_array[:N_atoms] + ipd_energies
        χ_AU_difference = χ_AU_array[i+1] - χ_AU_array[i]
        ionization_ratio = saha_equation(ne, T, χ_AU_difference, degeneracy_ratio) # n_{i+1}/n_{i}
        unnormalized_fractions.append(ionization_ratio*unnormalized_fractions[-1])

    xi_array = unnormalized_fractions/np.sum(unnormalized_fractions)
    return xi_array, χ_AU_array # ionization fractions and their energies
    
def calculate_ionization_fractions(plasma, nn, T, initial_guess = None, IPD = True, N_ions=None):
	# Initial guess for ionization fractions (all neutral initially)
	Z, χ0_AU_array, g_degeneracy_array = plasma.Z, plasma.χ0_AU_array, plasma.g_degeneracy_array
	if N_ions is None:
		N_atoms = len(χ0_AU_array)
		N_ions = N_atoms - 1
	else:
		N_atoms = N_ions + 1

	if initial_guess is None:
		More_Zbar = More_TF_Zbar(Z, nn, T)
		fake_ionization_fraction = np.zeros(N_atoms)
		fake_ionization_fraction[1] = More_Zbar
		initial_guess = get_ionization_fractions(plasma, fake_ionization_fraction, nn*More_Zbar, T, IPD=False)[0].copy()

	fractions = np.array(initial_guess).copy()
	# Minimize the residuals to find the self-consistent ionization fractions

	α = 0.1
	iters = 0
	while iters < 1000:
		Zbar_old = get_Zbar(fractions)
		new_fractions, χs = get_ionization_fractions(plasma, fractions, nn, T, IPD=IPD)
		fractions = α*new_fractions + (1-α)*fractions
		Zbar = get_Zbar(fractions)
		err = Zbar/Zbar_old-1
		# print(f"Zbar Old: {Zbar_old:0.3f}, New: {Zbar:0.3f}, rel change: {err:0.3e}")
		iters +=1
		# print("\tMinimizer success = ", result)
		# for i, fraction in enumerate(fractions):
		#     print(f"Ionization state {i}: {fraction:.2e}")
		    
	# Calculate the average ionization
	average_ionization = get_Zbar(fractions)
	return average_ionization, fractions, χs

