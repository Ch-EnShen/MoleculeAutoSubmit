"""
This program is to write Gaussian input file automatically.
- No solvent

Parameters
- molecule_name
- state
- job_type
- functional
- basis_set
- memeory
- core_num
"""

from pathlib import Path
import pickle
import os

# choices of functional and basis set
functional_arr = ["b3lyp", "pbe1pbe"]
basis_set_arr = ["6-31g(d)", "6-311+g(2d,2p)"]

# set required parameters
operating_sys = "windows"
molecule_name = "PTM"
state = "s0"
job_type = "opt-freq"
functional = functional_arr[0]
basis_set = basis_set_arr[1]
file_name = f"{molecule_name}-{functional}-{state}-{job_type}"
memory = 192 #GB
core_num = 32

# create .com file
Gaussian_input = Path(f"{molecule_name}/{file_name}.gjf")
Gaussian_input.touch(exist_ok = True)
f = open(Gaussian_input, 'w')

# write .com file
f.write(f"%chk={file_name}.chk\n")
f.write(f"%mem={memory}GB\n")
f.write(f"%nprocshared={core_num}\n")

# job type
if job_type == "opt-freq" and state == 's0':
    input_line = f"#p opt=(calcfc,tight) freq {functional}/{basis_set} int=ultrafine nosymm\n\n"
elif job_type == "opt-freq" and state == 's1':
    input_line = f"#p opt=(calcfc,tight) freq td {functional}/{basis_set} int=ultrafine nosymm\n\n"
elif job_type == "nacme":
    input_line = f"#p td {functional}/{basis_set} prop=(fitcharge,field) iop(6/22=-4, 6/29=1, 6/30= 0, 6/17=2) nosymm\n\n"
# elif job_type == "td":
    # input_line = f""

f.write(input_line)
# For windows: Chem3D
if operating_sys == "windows":
    f.write("Title\n\n")
# For Mac: Avagadro
if operating_sys == "Mac":
    f.write("Title\n")

# input geometry
## from .gjf
structure = open(f"{molecule_name}/{molecule_name}.gjf", 'r')
f.writelines(structure.readlines()[4:])

f.close()
