import pubchempy as pcp
from urllib.request import urlopen
from urllib.parse import quote
from rdkit import Chem
from rdkit.Chem import AllChem

def name2SMILES(compound_name: str) -> str:
    """ 
    Generte SMILES from compound name.
    The compound is searched in the NCI/CADD Chemical Identifier Resolver.
    (website: https://cactus.nci.nih.gov/chemical/structure)
    If error, use the PubChem PUB REST API to search the compound.
    (website: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
    (module: https://github.com/mcs07/PubChemPy)
    Fetching data using the urllib library.
    (tutorial: https://docs.python.org/3/howto/urllib2.html)
    
    Args:
        compound_name (str): Name of the compound
    Return: 
        str: SMILES string of the compound
    """
    # code ref: https://stackoverflow.com/questions/54930121/

    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/'\
            + quote(compound_name)\
            + '/smiles'
        SMILES = urlopen(url).read().decode('utf8')
        return SMILES
        
    except Exception as e:
        print(f"Compound \"{compound_name}\" is not found " +\
              f"in the NIH Chemical Identifier Resolver.\n" +\
              f"Error:\n{e}\n")
        try:
            possible_compounds = pcp.get_compounds(f"{compound_name}", 'name')
            if (len(possible_compounds) == 1):
                return possible_compounds[0].isomeric_smiles
            else:
                for i, pc in enumerate(possible_compounds):
                    print(
                        f"Possible Compound {i+1}: " +\
                        f"{pc.iupac_name}(CID: {pc.cid})\n")
                c_num = "Multiple" if len(possible_compounds) > 1 else "No"
                raise TypeError(
                    f"{c_num} compounds found for \"{compound_name}\".\n")
        except:
            raise ValueError(f"Compound \"{compound_name}\" is not found.\n")

def calculation_option_config(calculation_option:list):
    input_options = []
    for option in calculation_option:
        if option == "opt":
            input_options.append("opt=(calcfc,tight)")
        elif option == "freq":
            input_options.append("freq")
    return " ".join(input_options)

def gjf_name(
    mol_name:str, 
    functional: str,
    state:str, 
    calculation_options:list
) -> str:
    """ Generate Gaussian input file name """
    return f"{mol_name}-{functional}-{state}-{'-'.join(calculation_options)}.gjf"

def SMILES2gjf(
    SMILES: str, 
    mol_name: str,
    memory: int = 192,
    nprocshared: int = 32,
    state: str = "s0",
    calculation_options: list = ["opt", "freq"],
    functional: str = "b3lyp",
    basis: str = "6-31g++(d,p)"
) -> None:
    
    """ Generate gjf file from SMILES """
    mol_name = mol_name.replace(" ", "_")
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))

    if mol is not None:
        # Generate 3D coordinates for the molecule
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Write the Gaussian input file (.gjf)
        calc_option_config_str = calculation_option_config(calculation_options)

        output_file_name = gjf_name(
            mol_name, functional, state, calculation_options)
        with open(output_file_name, 'w') as f:
            f.write(f'%chk={output_file_name.replace(".gjf", ".chk")}\n')
            f.write(f'%mem={memory}GB\n')
            f.write(f'%nprocshared={nprocshared}\n')
            f.write(f'#p {calc_option_config_str} {functional}/{basis} int=ultrafine nosymm\n\n')
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

def name2gjf(compound_name:str):
    """ Generate gjf file from compound name """
    SMILES2gjf(
        SMILES=name2SMILES(compound_name), mol_name=compound_name)
    
    return None

if __name__ == "__main__":
    """ testing
    names  = ['3-Methylheptane', 'Aspirin', 'Diethylsulfate', 'Diethyl sulfate', '9-cyano-anthracene', 'Adamant']
    for name in names:
        name2gjf(name)
    """
    names = ['9-cyano-anthracene']
    for name in names:
        name2gjf(name)
    # names = sys.argv[1:]