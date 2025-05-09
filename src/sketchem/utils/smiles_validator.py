# This file contains the main function that compares a drawn molecule against a target SMILES string using AI


from rdkit import Chem
from rdkit.Chem import rdFMCS
from google import genai
from google.genai import types


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
