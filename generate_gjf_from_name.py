import sys
from urllib.request import urlopen
from urllib.parse import quote
from rdkit import Chem
from rdkit.Chem import AllChem

def Name2SMILES(compound_name:str):
    """ Generte SMILES from compound name """
    # ref: https://stackoverflow.com/questions/54930121/converting-molecule-name-to-smiles/54932071#54932071?newreg=6e3c8e50485b4857bc8639d31bfdd65c

    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(compound_name) + '/smiles'
        SMILES = urlopen(url).read().decode('utf8')
        return SMILES
        
    except:
        raise ValueError(f"Compound \"{compound_name}\" is not found.\n")

def generate_calcaulation_options_str(calculation_options:list):
    """ Generate calculation options string """
    calculation_options_str = []
    if "opt" in calculation_options:
        calculation_options_str.append("opt")
    if "freq" in calculation_options:
        calculation_options_str.append("freq")
    return "-".join(calculation_options)

def calculation_option_config(calculation_option:list):
    input_options = []
    for option in calculation_option:
        if option == "opt":
            input_options.append("opt=(calcfc,tight)")
        elif option == "freq":
            input_options.append("freq")
    return " ".join(input_options)

def gjf_name(mol_name:str, state:str, calculation_options:list):
    """ Generate Gaussian input file name """
    calculation_options_str = generate_calcaulation_options_str(calculation_options)
    return f"{mol_name}-{state}-{calculation_options_str}.gjf"

def SMILES2gjf(SMILES:str, 
               mol_name:str,
               memory = 192,
               nprocshared = 32,
               state = "s0",
               calculation_options = ["opt", "freq"],
               functional = "b3lyp",
               basis = "6-31g++(d,p)"):
    
    """ Generate gjf file from SMILES """
    mol_name = mol_name.replace(" ", "_")
    mol = Chem.MolFromSmiles(SMILES)

    if mol is not None:
        # Generate 3D coordinates for the molecule
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Write the Gaussian input file (.gjf)
        calc_option_config_str = calculation_option_config(calculation_options)

        output_file_name = gjf_name(mol_name, state, calculation_options)
        with open(output_file_name, 'w') as f:
            f.write(f'%chk={output_file_name.replace(".gjf", ".chk")}\n')
            f.write(f'%mem={memory}GB\n')
            f.write(f'%nprocshared={nprocshared}\n')
            f.write(f'# {calc_option_config_str} {functional}/{basis} int=ultrafine nosymm\n\n')
            f.write(Chem.MolToXYZBlock(mol))
            f.write("\n\n")

        """
        Update the title line (in line 6) and spin multiplicity (default is 0 1)
        in the Gaussian input file
        """
        ### start ###
        with open(output_file_name, 'r') as f:
            lines = f.readlines()
        new_content = f"{mol_name}/{state}\n\n" + "  0  1"
        lines[6 - 1] = new_content
        with open(output_file_name, 'w') as f:
            f.writelines(lines)
        ### end ###
        
        print(f"Conversion successful. Gaussian input file '{mol_name}' generated.")
    else:
        raise ValueError(f"Error for molecule {mol_name}: \
                           Unable to convert SMILES to molecule object in RDkit.")

def Name2gjf(compound_name:str):
    """ Generate gjf file from compound name """
    SMILES2gjf(SMILES   = Name2SMILES(compound_name), 
               mol_name = compound_name)
    
    return None

if __name__ == "__main__":
    """ testing
    names  = ['3-Methylheptane', 'Aspirin', 'Diethylsulfate', 'Diethyl sulfate', '9-cyano-anthracene', 'Adamant']
    for name in names:
        Name2gjf(name)
    """
    names = sys.argv[1:]
    for name in names:
        Name2gjf(name)