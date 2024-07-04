import os
import numpy as np

# Get the current folder path
current_folder = os.getcwd()

# Iterate through each subfolder in the current folder
for subfolder in os.listdir(current_folder):
    if os.path.isdir(subfolder):
        subfolder_path = os.path.join(current_folder, subfolder)
        if os.path.isfile(subfolder_path + "/adiabatic_energy.npy"):
            # Load the data from "data.npy"
            adiabatic_energy = np.load(subfolder_path + "/adiabatic_energy.npy")
            # Print the data
            print(f"{subfolder}:{adiabatic_energy} eV")
