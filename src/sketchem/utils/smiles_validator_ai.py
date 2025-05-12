# This file contains the main function that compares a drawn molecule against a target SMILES string using AI


import streamlit as st
from rdkit import Chem
from rdkit.Chem import rdFMCS
from google import genai
from google.genai import types
from sketchem.utils.smiles_validator_ai_prompt import ai_prompt
from streamlit.logger import get_logger
import logging

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)


def validate_drawing_with_ai(api_key, image_bytes: bytes, target_smiles: str, threshold: float = 0.85) -> bool | str:
    # Check for empty API key first
    if not api_key:
        return "❗ Gemini API key not set."
    
    try:

        # Call the Gemini API to get the smiles 
        client = genai.Client(api_key=api_key)
        prompt = ai_prompt
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[
                types.Part.from_bytes(data=image_bytes, mime_type="image/png"),
                prompt,
            ],
        )
        logger.info(f"Gemini Detected Mol: {response.text.strip()}")
        st.session_state.last_gemini_detected_mol = response.text.strip()
        
        if response.text.strip() == "INVALID_STRUCTURE":
            return response.text.strip()
        
        # Re-use MCS code from decimer integration
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
