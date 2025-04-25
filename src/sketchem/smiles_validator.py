# This file contains the main function that compares a drawn molecule against a target SMILES string using various comparison methods (the only usefuls one are morgan and mcs, the rest are POCs)

import DECIMER
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import numpy as np
from rdkit.Chem import rdFMCS

def validate_drawing(image_path: str, target_smiles: str, method: str = 'mcs', threshold: float = 0.85) -> bool:
    """
    Validates a drawn molecule against a target SMILES string using various comparison methods.
    
    Args:
        image_path (str): Path to the image file containing the drawn molecule
        target_smiles (str): SMILES string to compare against
        method (str): Comparison method ('morgan', 'exact', 'canonical', or 'mcs')
        threshold (float): Similarity threshold for Morgan fingerprint or MCS comparison
    
    Returns:
        bool: True if the molecules match according to the specified method
    """
    try:
        # Get SMILES from image using DECIMER
        model_name = "Canonical"

        # Returned as a tuple that looks like ('CC[O]', [('C', 0.7732562), ('C', 0.8388427), ('[', 0.5454414), ('O', 0.9842202), (']', 0.6354147)]) hence need only the first element
        predicted_smiles = DECIMER.predict_SMILES(image_path, model_name, hand_drawn = True)[0] 
        print(f"Predicted SMILES: {predicted_smiles}") 


        # Method 1: Morgan Fingerprint Similarity -> broken
        if method == 'morgan':
            # Generate molecules
            mol1 = Chem.MolFromSmiles(predicted_smiles)
            mol2 = Chem.MolFromSmiles(target_smiles)
            
            # Check whether the molecules exist
            if mol1 is None or mol2 is None:
                return False

            # Generate Morgan fingerprints
            fpgen = GetMorganGenerator(radius=2, fpSize=2048)
            fp1 = fpgen.GetFingerprint(mol1)
            fp2 = fpgen.GetFingerprint(mol2)
            
            # Calculate Tanimoto similarity
            # Threshold inputted only works here
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            return similarity >= threshold 
            
        
        # Method 2: Exact SMILES match (not exactly recommended)
        elif method == 'exact':
            return predicted_smiles == target_smiles
        
        # Method 3: Canonical SMILES comparison
        elif method == 'canonical':
            # Convert both to molecules 
            mol1 = Chem.MolFromSmiles(predicted_smiles)
            mol2 = Chem.MolFromSmiles(target_smiles)
            
            if mol1 is None or mol2 is None:
                return False
            
            #Convert back to canonical SMILES  but canonical
            canon1 = Chem.MolToSmiles(mol1, canonical=True)
            canon2 = Chem.MolToSmiles(mol2, canonical=True)
            
            return canon1 == canon2
            
        # Method 4: Maximum Common Substructure comparison
        elif method == 'mcs':
            mol1 = Chem.MolFromSmiles(predicted_smiles)
            mol2 = Chem.MolFromSmiles(target_smiles)
            
            if mol1 is None or mol2 is None:
                return False
            
            # Find MCS
            mcs = rdFMCS.FindMCS([mol1, mol2], 
                                atomCompare=rdFMCS.AtomCompare.CompareElements,
                                bondCompare=rdFMCS.BondCompare.CompareOrder,
                                matchValences=True,
                                ringMatchesRingOnly=True,
                                completeRingsOnly=True)
            
            # Calculate similarity based on MCS size
            if mcs.numAtoms == 0:
                return False
                
            # Calculate similarity as fraction of atoms in common
            similarity = mcs.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
            return similarity >= threshold
            
        else:
            raise ValueError("Method must be 'morgan', 'exact', 'canonical', or 'mcs'")
            
    except Exception as e:
        print(f"Error inside function validate_drawing: {str(e)}")
        return False

"""
Recommended thresholds:
For Morgan method:
-need to run tests

For MCS method:
-need to run tests
"""
