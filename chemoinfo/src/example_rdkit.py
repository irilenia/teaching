from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
import numpy as np

results={"lipinski_true":0,
        "lipinski_false":0}

molecules = Chem.SDMolSupplier('data/PubChem_substance_text_covid-19_records.sdf')
print("Number of molecules read from file= ", len(molecules))

#this chunk goes through the molecules read in from the sdf file
#and calculates descriptors linked to the Lipinski rule of 5s. It
#then prints out which molecules are not Lipinski-compliant.
for mol in molecules:
    molecular_weight = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_bond_donor = Descriptors.NumHDonors(mol)
    h_bond_acceptors = Descriptors.NumHAcceptors(mol)
    if molecular_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 5 :
        results["lipinski_true"] += 1
    else:
        results["lipinski_false"] += 1
        print("Not Lipinski compliant: ", mol.GetProp('_Name'))
print(results)


#this chunk fingerprints the molecules (and prints out one of the fingerprints as an example
#and then calculates their similarity,
#printing out a matrix of Tanimoto scores for all comparisons.
list_of_fps = [Chem.RDKFingerprint(i) for i in molecules]
matrix = np.zeros([len(list_of_fps),len(list_of_fps)])

for i in range(len(list_of_fps)):
    for j in range(len(list_of_fps)):
        tani = DataStructs.FingerprintSimilarity(list_of_fps[i],list_of_fps[j])
        matrix[i][j] = tani
print(list_of_fps[0].ToBitString())
print(matrix)
