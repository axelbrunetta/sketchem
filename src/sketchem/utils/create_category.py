# This file contains the main function that creates a new molecule category using Gemini AI and smiles from Pubchem 

from google import genai
import pubchempy as pcp
from streamlit.logger import get_logger
import logging
import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from sketchem.utils.create_category_prompt import create_prompt

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def check_category_is_default(selected_category):
    if selected_category in MOLECULE_CATEGORIES.keys():
        st.session_state.category_is_default = True
    else:
        st.session_state.category_is_default = False

def get_molecules_for_category_pubchem(api_key, user_prompt, jupyternb: bool = False):
    full_prompt = create_prompt(user_prompt)

    try:
        # Check for empty API key first
        if not api_key:
                return "Gemini API key not set."
        
        # Call the Gemini API to get the category
        client = genai.Client(api_key=api_key)
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[full_prompt],
        )
        response_text = response.text.strip()
        
        lines = response_text.split('\n')
        if not lines:
            return "Error: AI returned empty response.", []

        # First line is expected to be "Category Name (number)"
        category_header = lines[0]
        molecule_names_from_ai = [name.strip() for name in lines[1:] if name.strip()]

        if not molecule_names_from_ai:
            return f"Error: AI returned category '{category_header}' but no molecule names.", []

        logger.info(f"Created category (names only): {category_header}, Elements: {molecule_names_from_ai}")

        molecules_with_smiles = []
        for name in molecule_names_from_ai:
            try:
                compounds = pcp.get_compounds(name, 'name')
                if compounds:
                    # Take the first compound found
                    compound = compounds[0]
                    if compound.canonical_smiles:
                        molecules_with_smiles.append({
                            "name": name, # Or compound.iupac_name if preferred
                            "smiles": compound.canonical_smiles
                        })
                        if not jupyternb:
                            logger.info(f"Found SMILES for '{name}': {compound.canonical_smiles}")
                    else:
                        if not jupyternb:
                            logger.info(f"Warning: No canonical SMILES found for '{name}' in PubChem.")
                        # just continue
                        continue
                else:
                    if not jupyternb:
                        logger.info(f"Warning: No molecule found for '{name}' in PubChem.")
                    # skip this molecule
                    continue
            except pcp.PubChemPyError as e:
                if not jupyternb:
                    logger.info(f"PubChemPy Error for '{name}': {e}")
                # Skip this molecule 
                continue
            except Exception as e:
                if not jupyternb:
                    logger.info(f"General Error processing '{name}': {e}")
                # Skip this molecule
                continue
        
        # Add the new category to the additional_categories state var
        molecules_dict = {item["name"]: item["smiles"] for item in molecules_with_smiles}

        if not jupyternb:
            st.session_state.additional_categories[category_header] = molecules_dict
            
            # Store the newly created category name temporarily to select it later
            st.session_state.last_created_category = category_header
            
            # Immediately update the selected category and set category_is_default to False
            st.session_state.selected_molecule_category = category_header
            st.session_state.category_is_default = False
            
            logger.info(f"Created category: {category_header}, Elements: {molecules_dict}")
            logger.info(f"Updated category: {st.session_state.selected_molecule_category}, isDefault: {st.session_state.category_is_default}")
        
            return "Successfully created category"
        else:
            return category_header, molecules_dict

    except Exception as e:
        logger.info(f"Error communicating with Gemini or processing response: {e}")
        return f"Error: {e}", []

