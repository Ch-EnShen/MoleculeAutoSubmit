# MoleculeAutoSubmit

Submit Gaussian job on PBS (Portable Batch System) clusters by giving only molecule names.

## Usage
Run the following command on a PBS environment
```
./functional_submission.sh
```
Enter the following adjustments:
```
Queue? (Defaults to mem192v2)
No. of Nodes? (Defaults to 1)
Jobname? (Defaults to mol_sub_default)
Functional? (Defaults to b3lyp)
```

## Seperate Scripts
### generate_gjf_from_chemical_name
generate Gaussian input file (.gjf) from chemical name
#### Usage: 
```
python generate_gjf_from_name.py [MOLECULAER_NAMES]
```
for example, 
```
python generate_gjf_from_name.py aspirin anthracene
```

### runscript_get_dipole_and_energy
parse dipole moment and adiabatic energy from log files of Gaussian calculation with `opt freq` keywords.
#### Usage:
```
python runscript_get_dipole_and_energy.py [MOLECULAER_NAME] [FUNCTIONAL]
```

### runscript_mol_sub
generate s1-opt-freq.gjf file from s0-opt-freq.log and .gjf files
generate s1-nacme.gjf file from s1-opt-freq.log and .gjf files
#### Usage
```
python runscript_mol_sub.py [MOLECULAER_NAME] [CALC_TYPE] [FUNCTIONAL]
```
with `[CALC_TYPE] = opt-freq`, the s1-opt-freq.gjf file is generated.
with `[CALC_TYPE] = nacme`, the s1-nacme.gjf file is generated.
