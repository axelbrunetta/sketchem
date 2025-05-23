"""
AI-powered molecule validation for Sketchem.

This file contains functions that use Google's Gemini AI to analyze
hand-drawn molecular structures and validate them against target molecules.
"""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import rdFMCS
from google import genai
from google.genai import types
import pubchempy as pcp
from sketchem.utils.smiles_validator_ai_prompt import ai_prompt
from streamlit.logger import get_logger
import logging

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def validate_drawing_with_ai(response, target_smiles, threshold, jupyternb: bool = False):
    """
    Validate Gemini's response against the target molecule.
    
    Args:
        response: The Gemini API response containing the molecule names
        target_smiles: The target molecule's SMILES string
        threshold: The MCS similarity threshold for validation
        jupyternb: Whether we are running in a Jupyter notebook (for testing)
        
    Returns:
        bool: True if the drawing is valid, False otherwise
    """
    # Re-use MCS code from decimer integration
    try:
        # Parse the response to get the three names
        names = response.text.strip().split('\n')
        
        # Filter out empty names
        names = [name for name in names if name.strip()]
        
        if not names:
            logger.error("No valid names returned from Gemini")
            return False
            
        # For each name, try to get SMILES from PubChem
        for name in names:
            try:
                compounds = pcp.get_compounds(name, 'name')
                if not compounds:
                    if not jupyternb:
                        logger.info(f"No compounds found for name: {name}")
                    continue
                if not jupyternb:
                    logger.info(f"Found compound SMILES for {name}: {compounds[0].canonical_smiles}")
                # Get the first compound's SMILES
                mol1 = Chem.MolFromSmiles(compounds[0].canonical_smiles)
                mol2 = Chem.MolFromSmiles(target_smiles)
                
                if mol1 is None or mol2 is None:
                    logger.info(f"Invalid SMILES for name: {name}")
                    continue
                
                # MCS-based similarity 
                mcs = rdFMCS.FindMCS([mol1, mol2], 
                                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                                    bondCompare=rdFMCS.BondCompare.CompareOrder,
                                    matchValences=True,
                                    ringMatchesRingOnly=True,
                                    completeRingsOnly=True)
                
                if mcs.numAtoms == 0:
                    continue
                
                mcs_similarity = mcs.numAtoms / max(mol1.GetNumAtoms(), mol2.GetNumAtoms())
                
                # If this name's molecule matches above threshold, return True
                if mcs_similarity >= threshold:
                    if not jupyternb:
                        logger.info(f"Match found for name: {name} with similarity: {mcs_similarity}")
                    return True
                    
            except Exception as e:
                logger.error(f"Error processing name {name}: {e}")
                continue
        
        # If we get here, none of the names matched
        return False
    except Exception as e:
        logger.error(f"Error in validate_drawing_with_ai: {e}")
        return False

def get_molecule_with_ai(api_key, image_bytes: bytes, target_smiles: str, threshold: float = 0.85, jupyternb: bool = False) -> bool | str:
    """
    Validate a hand-drawn molecule against a target SMILES string using Gemini AI.
    
    Args:
        api_key: The Gemini API key
        image_bytes: The image bytes containing the hand-drawing
        target_smiles: The target molecule's SMILES string
        threshold: The MCS similarity threshold for validation
        jupyternb: Whether we are running in a Jupyter notebook (for testing)
        
    Returns:
        bool: True if the drawing is valid, False otherwise
        str: Error message if there was an error
    """
    # Check for empty API key first
    if not api_key:
        return "❗ Gemini API key not set."
    
    try:
        # Call the Gemini API to get the molecule names
        client = genai.Client(api_key=api_key)
        prompt = ai_prompt
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[
                types.Part.from_bytes(data=image_bytes, mime_type="image/png"),
                prompt,
            ],
        )
        
        st.session_state.last_gemini_detected_mol = response.text.strip()
        
        if not jupyternb:
            logger.info(f"Gemini Detected Molecule Names: {response.text.strip()}")
            return validate_drawing_with_ai(response, target_smiles, threshold)
        else:
            return response
    except Exception as e:
        return f"❗ Gemini API error: {e}"
