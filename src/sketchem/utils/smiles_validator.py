# This file contains the main function that compares a drawn molecule against a target SMILES string using various comparison methods (the only usefuls one are mcs and using gemini, the rest are POCs)

import DECIMER
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem import rdFMCS
from google import genai
from google.genai import types
from PIL import Image
import io


def validate_drawing(image_path: str, target_smiles: str, method: str = 'mcs', threshold: float = 0.85) -> bool:
    """
    Validates a drawn molecule against a target SMILES string using various comparison methods.
    
    Arguments:
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


        # Method 1: Morgan Fingerprint Similarity -> broken since decimer returns things like CCO.[I-] for ethanol -> go with method 4 or gemini method
        """
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
         """   
        
        # Method 2: Exact SMILES match (not exactly recommended) -> Proof of concept, not useful
        """  
        elif method == 'exact':
            return predicted_smiles == target_smiles
        """  
        # Method 3: Canonical SMILES comparison -> Proof of concept, not useful
        """  
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
        """   
        # Method 4: Maximum Common Substructure comparison -> works even when decimer makes errors of interpretation such as CCO.[I-] for ethanol
        if method == 'mcs':
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


def validate_drawing_with_ai(api_key, image_bytes: bytes, target_smiles: str, threshold: float = 0.85) -> bool:
    # Check for empty API key first
    if not api_key:
        return "❗ Gemini API key not set."
    
    try:

        # Call the Gemini API to get the smiles 
        client = genai.Client(api_key=api_key)
        prompt = """

You are an expert chemist and molecular recognition system. Your task is to carefully analyze hand-drawn chemical structures and provide the correct SMILES (Simplified Molecular Input Line Entry System) notation for the molecule.

It is extremely important that you get the SMILES string exactly right. Double-check for correctness.

Follow these guidelines:
	1.	Identify all atoms and bonds precisely from the drawing, including:
	•	Single, double, triple bonds
	•	Aromatic rings (benzene-like structures)
	•	Branches and substituents
	2.	Do not assume or guess any part of the structure. Only provide the SMILES for what is clearly depicted.
	3.	Verify your SMILES string by mentally reconstructing the molecule to ensure it matches the drawing.

⸻

Examples:
	1.	Drawing: A benzene ring (hexagon with alternating double bonds)
Correct SMILES: c1ccccc1
	2.	Drawing: Ethanol (CH3-CH2-OH)
Correct SMILES: CCO
	3.	Drawing: Cyclohexane (hexagon with single bonds only)
Correct SMILES: C1CCCCC1
	4.	Drawing: Acetic acid (CH3-COOH)
Correct SMILES: CC(=O)O

⸻

IMPORTANT: After providing the SMILES string, do not include any other explanations or commentary. Simply output the SMILES code.
If it is not a valid molecule drawing, return "Invalid molecule drawing"

"""
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[
                types.Part.from_bytes(data=image_bytes, mime_type="image/png"),
                prompt,
            ],
        )
        if response.text.strip() == "Invalid molecule drawing":
            return False
        
        # Re-use MCS code from above
        mol1 = Chem.MolFromSmiles(response.text.strip())
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
        
        if mcs.numAtoms == 0:
            return False
            
        # Calculate similarity as fraction of atoms in common
        similarity = mcs.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
        return similarity >= threshold
    
    except Exception as e:
        # This should throw type error since wrong type
        return f"❗ Gemini API error: {e}"
