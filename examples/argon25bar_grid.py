# Zach Johnson 6/27/24
# Xenon 5 bar experiment from UCLA

import numpy as np
from hnc.hnc.constants import *

from saha.core.saha import plasma
from saha.core.table_generator import saha_table

nn_invcc_at_Pbar_TK = lambda Pbar, TK: Pbar*bar_to_AU/(TK*K_to_AU)*AU_to_invcc

Ar_nn_invcc = nn_invcc_at_Pbar_TK(25, 290)
Ar_TK_peak = 17.761029411764707e3 # 0.008097165991902834 ns?

# https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=Ar&units=1&at_num_out=on&el_name_out=on&seq_out=on&shells_out=on&level_out=on&e_out=0&unc_out=on&biblio=on
Ar_ionization_energies_AU = np.array([0, 15.7596119, 27.62967, 40.735, 59.58, 74.84, 91.290])*eV_to_AU
Ar_J_ground_level = np.array([0, 3/2, 2, 3/2, 0, 1/2, 0])
Ar_ionization_degeneracies = 2*Ar_J_ground_level + 1

# Input: Number of ionizations
Z = 18
Ar25bar_plasma = plasma("Ar25bar", Z, Ar_ionization_energies_AU, Ar_ionization_degeneracies)

# Temperature density
nn_AU_range = Ar_nn_invcc*invcc_to_AU * np.array([0.1, 5])  # 1/cc
T_AU_range = Ar_TK_peak * K_to_AU * np.array([0.1, 5])# Kelvin

# Calculate ionization fractions

Ar25bar_table = saha_table(Ar25bar_plasma, nn_AU_range, T_AU_range, Nn=20, NT=200, N_ions = 6)

