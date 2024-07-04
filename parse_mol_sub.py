import numpy as np
import os
import copy
from pathlib import Path
import pickle

def parse_log(filename):
    with open(filename,'r') as f:
        total_file = f.readlines()
        total_file = reversed(total_file)
        
        start = 0
        geometry_coordinate = []
        for line in total_file:
            line = line.strip()
            if "Distance matrix (angstroms):" in line:
                start += 1
                continue
            elif start > 0 and "---------------------------" in line:
                start += 1
                continue                
            if start == 2:
                geometry_coordinate.append(line.split()[3:])
            elif start > 2:
                break
            
    return geometry_coordinate[::-1]

def get_trans_dip_adiabatic_energy(s1_log, s0_log):
    """
    Args:
        s1_log: <molecule>-opt-freq-s1.log
        s0_log: <molecule>-opt-freq-s0.log
    Return:
        transition dipole vector (Debye)
        transition dipole magnitude (Debye)
        Adiabatic energy (eV)
    """
    with open(s1_log, 'r') as f:
        s1_log_file = f.readlines()
        reverse_s1_log_file = reversed(s1_log_file)
        count = 0
        
        lines = []
        for line in reverse_s1_log_file:
            if len(lines) < 3:
                lines.append(line)
                continue
            else:
                lines.pop(0)
                lines.append(line)
            
            if "Ground to excited state transition electric dipole" in lines[-1]:
                x, y, z, dip_square = (lines[0].split())[1:5]
                trans_dip_vector = np.array([float(x) , float(y), float(z)]) * 2.5416
                trans_dip_magnitude = np.sqrt(float(dip_square)) * 2.5416
                break
        
        for line in reverse_s1_log_file:
            if "Total Energy" in line:
                s1_minimun_energy_Hartree = float(line.split()[-1])
                break
    
    with open(s0_log, 'r') as f:
        s0_log_file = f.readlines()
        for line in reversed(s0_log_file):
            if "SCF Done" in line:
                s0_minimun_energy_Hartree = float(line.split()[4]) 
                break
    
    adiabatic_energy = abs(s1_minimun_energy_Hartree - s0_minimun_energy_Hartree) * 27.2107
    return trans_dip_vector, trans_dip_magnitude, adiabatic_energy

def write_in_gjf_s1(log_filename:str, s0_gjf_filename:str, molecule_name:str):
    current_dir = os.getcwd()
    s1_init_geometry_coordinate = parse_log(current_dir + f"/{molecule_name}/{log_filename}")
    # print(s1_init_geometry_coordinate)
    # touch s1.gjf file
    
    s1_gjf_filename = s0_gjf_filename.replace("s0", "s1")
    s1_gjf = Path(current_dir + f"/{molecule_name}/{s1_gjf_filename}")
    s1_gjf.touch(exist_ok = True)

    s0_gjf = open(current_dir + f"/{molecule_name}/{s0_gjf_filename}", 'r')
    s0_gjf_list = s0_gjf.readlines()
    # print(s0_gjf_list)
    with open(s1_gjf, 'w') as f:
        for i in range(8):
            if i == 0:
                f.write(s0_gjf_list[i].replace("s0", "s1"))
            elif i == 3:
                job_input_line = s0_gjf_list[3].split("freq")#.insert(1, "freq td")
                job_input_line.insert(1, "freq td")
                f.write("".join(job_input_line))
            else:
                f.write(s0_gjf_list[i])
        j = 8
        try:
            while True:
                line = s0_gjf_list[j]
                newline = copy.deepcopy(line)
                start_splitting = 2 if target_geometry_coordinate[j-8][1] == '0' else 1
                for i, num in enumerate(line.split()[start_splitting:]):
                    newline = newline.replace(num, s1_init_geometry_coordinate[j-8][i])
                    # print(newline)
                f.write(newline)
                j += 1
                # print(j)
        except IndexError as e:
            print(e)
    s0_gjf.close()
    return 

def write_in_gjf(input_log_filename:str, input_gjf_filename:str, molecule_name:str, functional:str, calc_type:str):
    current_dir = os.getcwd()
    target_geometry_coordinate = parse_log(current_dir + f"/{molecule_name}_{functional}/{input_log_filename}")
    
    output_gjf_filename = input_gjf_filename.replace("s0", "s1") if calc_type == "opt-freq" else input_gjf_filename.replace("opt-freq", "nacme").replace("s0", "s1")
    output_gjf = Path(current_dir + f"/{molecule_name}_{functional}/{output_gjf_filename}")
    output_gjf.touch(exist_ok = True)
    
    input_gjf = open(current_dir + f"/{molecule_name}_{functional}/{input_gjf_filename}", 'r')
    input_gjf_list = input_gjf.readlines()
    with open(output_gjf, 'w') as f:
        for i in range(8):
            if i == 0:
                if calc_type == "opt-freq":
                    f.write(input_gjf_list[0].replace("s0", "s1"))
                else: #calc_type == nacme
                    f.write(input_gjf_list[0].replace("opt-freq", "nacme").replace("s0", "s1"))
            elif i == 3:
                if calc_type == "opt-freq":
                    job_input_line = input_gjf_list[3].split("freq")
                    job_input_line.insert(1, "freq td")
                    f.write("".join(job_input_line))
                elif calc_type == "nacme":
                    matches = [x for x in input_gjf_list[i].split() if '/' in x]
                    f.write(f"#p td {matches[0]} prop=(fitcharge,field) iop(6/22=-4, 6/29=1, 6/30= 0, 6/17=2) nosymm\n")
            else:
                f.write(input_gjf_list[i])
        
        j = 8
        try:
            while True:
                line = input_gjf_list[j]
                newline = copy.deepcopy(line)
                start_splitting = 2 if target_geometry_coordinate[j-8][1] == '0' else 1
                for i, num in enumerate(line.split()[start_splitting:]):
                    newline = newline.replace(num, target_geometry_coordinate[j-8][i])
                f.write(newline)
                j += 1
        except IndexError:
            pass
        f.write("\n\n")
    input_gjf.close()
    return 


if __name__ == "__main__":
    # data = parse_log("test.log")
    # print(data.shape)
    # write_in_gjf_s1("benzonitrile-b3lyp-s0-opt-freq.log", "benzonitrile-b3lyp-s0-opt-freq.gjf", "benzonitrile")
    # write_in_gjf("benzonitrile-b3lyp-s0-opt-freq.log", "benzonitrile-b3lyp-s0-opt-freq.gjf", "benzonitrile", "opt-freq")
    # write_in_gjf("benzonitrile-b3lyp-s1-opt-freq.log", "benzonitrile-b3lyp-s0-opt-freq.gjf", "benzonitrile", "nacme")
    trans_dip_vector, trans_dip_magnitude, adiabatic_energy = get_trans_dip_adiabatic_energy("penta_thiophene-b3lyp-s1-opt-freq.log", "penta_thiophene-b3lyp-s0-opt-freq.log")
    np.save("transition_dipole_orientation", trans_dip_vector/np.linalg.norm(trans_dip_vector))
    np.save("transition_dipole_magnitude", trans_dip_magnitude)
    np.save("adiabatic_energy", adiabatic_energy)

