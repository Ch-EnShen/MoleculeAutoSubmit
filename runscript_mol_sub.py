import sys
import parse_mol_sub

molecule_name = sys.argv[1]
calc_type = sys.argv[2]
funtional = sys.argv[3]

if calc_type == "opt-freq":
    s0_opt_log = molecule_name + f'-{funtional}-s0-opt-freq.log'
    s0_gjf = molecule_name + f'-{funtional}-s0-opt-freq.gjf'
    parse_mol_sub.write_in_gjf(s0_opt_log, s0_gjf, molecule_name, calc_type)
elif calc_type == "nacme":
    s1_opt_log = molecule_name + f'-{funtional}-s1-opt-freq.log'
    s1_gjf = molecule_name + f'-{funtional}-s1-opt-freq.gjf'
    parse_mol_sub.write_in_gjf(s1_opt_log, s1_gjf, molecule_name, calc_type)
else:
    print("Error calc_type")