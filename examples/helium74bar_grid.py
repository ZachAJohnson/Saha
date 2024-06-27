# Zach Johnson 6/27/24
# Helum 74 bar experiment from UCLA

import numpy as np
from hnc.hnc.constants import *

from saha.core.saha import plasma
from saha.core.table_generator import saha_table

nn_invcc_at_Pbar_TK = lambda Pbar, TK: Pbar*bar_to_AU/(TK*K_to_AU)*AU_to_invcc

He_nn_invcc = nn_invcc_at_Pbar_TK(74, 290)
He_TK_peak = 14.790528233151186e3 # 0.0031746031746031746 ns

# https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=He&units=1&at_num_out=on&el_name_out=on&seq_out=on&shells_out=on&level_out=on&e_out=0&unc_out=on&biblio=on
He_ionization_energies_AU = np.array([0, 24.587389011, 54.4177655282])*eV_to_AU
He_J_ground_level = np.array([0, 1/2, 0 ])
He_ionization_degeneracies = 2*He_J_ground_level + 1

# Input: Number of ionizations
Z = 2
He74bar_plasma = plasma("He74bar", Z, He_ionization_energies_AU, He_ionization_degeneracies)

# Temperature density
nn_AU_range = He_nn_invcc*invcc_to_AU * np.array([0.1, 5])  # 1/cc
T_AU_range = He_TK_peak * K_to_AU * np.array([0.1, 5])# Kelvin

# Calculate ionization fractions

He74bar_table = saha_table(He74bar_plasma, nn_AU_range, T_AU_range, Nn=20, NT=200, N_ions = None)

