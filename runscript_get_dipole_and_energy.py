import numpy as np
import parse_mol_sub
import sys

molecule_name = sys.argv[1]
functional = sys.argv[2]
s1_path = f"{molecule_name}_{functional}/{molecule_name}-{functional}-s1-opt-freq.log"
s0_path = f"{molecule_name}_{functional}/{molecule_name}-{functional}-s0-opt-freq.log"
trans_dip_vector, trans_dip_magnitude, adiabatic_energy = parse_mol_sub.get_trans_dip_adiabatic_energy(s1_path, s0_path)
np.save(f"{molecule_name}_{functional}/trans_dipole_orientation.npy", trans_dip_vector/np.linalg.norm(trans_dip_vector))
np.save(f"{molecule_name}_{functional}/trans_dipole_magnitude.npy", trans_dip_magnitude)
np.save(f"{molecule_name}_{functional}/adiabatic_energy.npy", adiabatic_energy)