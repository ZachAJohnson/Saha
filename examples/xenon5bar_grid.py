# Zach Johnson 6/27/24
# Xenon 5 bar experiment from UCLA

import numpy as np
from hnc.hnc.constants import *

from saha.core.saha import plasma
from saha.core.table_generator import saha_table

nn_invcc_at_Pbar_TK = lambda Pbar, TK: Pbar*bar_to_AU/(TK*K_to_AU)*AU_to_invcc

Xe_nn_invcc = nn_invcc_at_Pbar_TK(5, 290)
Xe_TK_peak = 16.60952380952381e3 # 0 ns?

# https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=Ar&units=1&at_num_out=on&el_name_out=on&seq_out=on&shells_out=on&level_out=on&e_out=0&unc_out=on&biblio=on
Xe_ionization_energies_AU = np.array([0, 12.1298437, 20.975, 31.05, 42.20 , 54.1, 66.703, 91.6, 105.9778, 179.84, 202.0, 229.02])*eV_to_AU
Xe_J_ground_level = np.array([0, 3/2, 2, 3/2, 0, 1/2, 0, 1/2, 0, 5/2, 4, 9/2])
Xe_ionization_degeneracies = 2*Xe_J_ground_level + 1

# Input: Number of ionizations
Z = 54
Xe5bar_plasma = plasma("Xe5bar", Z, Xe_ionization_energies_AU, Xe_ionization_degeneracies)

# Temperature density
nn_AU_range = Xe_nn_invcc*invcc_to_AU * np.array([0.9, 1.1])  # 1/cc
T_AU_range = Xe_TK_peak * K_to_AU * np.array([0.1, 5])# Kelvin

# Calculate ionization fractions

Xe5bar_table = saha_table(Xe5bar_plasma, nn_AU_range, T_AU_range, Nn=5, NT=20, N_ions = 2)

