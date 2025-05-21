#main function that compares drawn molecule with target SMILES string using gemini


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
    #reuse MCS code from decimer integration
    try:
        #parse response to get the three names
        names = response.text.strip().split('\n')
        
        #filter out empty names
        names = [name for name in names if name.strip()]
        
        if not names:
            logger.error("No valid names returned from Gemini")
            return False
            
        #for each name, try to get SMILES from PubChem
        for name in names:
            try:
                compounds = pcp.get_compounds(name, 'name')
                if not compounds:
                    if not jupyternb:
                        logger.info(f"No compounds found for name: {name}")
                    continue
                if not jupyternb:
                    logger.info(f"Found compound SMILES for {name}: {compounds[0].canonical_smiles}")
                #get SMILES of first compound
                mol1 = Chem.MolFromSmiles(compounds[0].canonical_smiles)
                mol2 = Chem.MolFromSmiles(target_smiles)
                
                if mol1 is None or mol2 is None:
                    logger.info(f"Invalid SMILES for name: {name}")
                    continue
                
                #MCS-based similarity 
                mcs = rdFMCS.FindMCS([mol1, mol2], 
                                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                                    bondCompare=rdFMCS.BondCompare.CompareOrder,
                                    matchValences=True,
                                    ringMatchesRingOnly=True,
                                    completeRingsOnly=True)
                
                if mcs.numAtoms == 0:
                    continue
                
                mcs_similarity = mcs.numAtoms / max(mol1.GetNumAtoms(), mol2.GetNumAtoms())
                
                #if molecule of name matches above threshold return True
                if mcs_similarity >= threshold:
                    if not jupyternb:
                        logger.info(f"Match found for name: {name} with similarity: {mcs_similarity}")
                    return True
                    
            except Exception as e:
                logger.error(f"Error processing name {name}: {e}")
                continue
        
        #if we get here: none of the names matched
        return False
    except Exception as e:
        logger.error(f"Error in validate_drawing_with_ai: {e}")
        return False

def get_molecule_with_ai(api_key, image_bytes: bytes, target_smiles: str, threshold: float = 0.85, jupyternb: bool = False) -> bool | str:
    #check for empty api key first
    if not api_key:
        return "❗ Gemini API key not set."
    
    try:
        #call gemini api to get molecule names
        client = genai.Client(api_key=api_key)
        prompt = ai_prompt
        
        #call gemini api
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[
                types.Part.from_bytes(data=image_bytes, mime_type="image/png"),
                prompt,
            ],
        )
        
        st.session_state.last_gemini_detected_mol = response.text.strip()
        
        #validate drawing
        if not jupyternb:
            logger.info(f"Gemini Detected Molecule Names: {response.text.strip()}")
            return validate_drawing_with_ai(response, target_smiles, threshold)
        else:
            return response
    except Exception as e:
        return f"❗ Gemini API error: {e}"