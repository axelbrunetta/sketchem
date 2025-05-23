"""
Category creation utilities for Sketchem.

This file contains functions to create custom molecule categories
using Gemini AI and PubChem data.
"""

from google import genai
import pubchempy as pcp
from streamlit.logger import get_logger
import logging
import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from sketchem.utils.create_category_prompt import create_prompt

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

#check if category is default
def check_category_is_default(selected_category):
    """
    Check if a category is one of the default categories.
    
    Updates the session state to indicate whether the selected category
    is a default category or a custom one.
    
    Args:
        selected_category: The name of the category to check
    """
    if selected_category in MOLECULE_CATEGORIES.keys():
        st.session_state.category_is_default = True
    else:
        st.session_state.category_is_default = False

#get molecules for category using pubchem
def get_molecules_for_category_pubchem(api_key, user_prompt, jupyternb: bool = False):
    """
    Generate a new molecule category based on user input.
    
    Uses Gemini AI to interpret the user's request and find relevant
    molecules from PubChem.
    
    Args:
        api_key: Google Gemini API key
        user_prompt: User's description of the desired category
        jupyternb: Whether this is being run in a Jupyter notebook
        
    Returns:
        Success message or error message
    """
    full_prompt = create_prompt(user_prompt)

    try:
        #check for empty api key first
        if not api_key:
                return "Gemini API key not set."
        
        #call Gemini api to get category
        client = genai.Client(api_key=api_key)
        
        #generate content
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[full_prompt],
        )
        response_text = response.text.strip()
        
        #split into lines
        lines = response_text.split('\n')
        if not lines:
            return "Error: AI returned empty response.", []

        #first line is expected to be "Category Name (number)"
        category_header = lines[0]
        molecule_names_from_ai = [name.strip() for name in lines[1:] if name.strip()]

        #error if no molecule names
        if not molecule_names_from_ai:
            return f"Error: AI returned category '{category_header}' but no molecule names.", []

        logger.info(f"Created category (names only): {category_header}, Elements: {molecule_names_from_ai}")

        #get smiles for molecules
        molecules_with_smiles = []
        for name in molecule_names_from_ai:
            try:
                compounds = pcp.get_compounds(name, 'name')
                if compounds:
                    #take the first compound found
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
            #error if pubchem error
            except pcp.PubChemPyError as e:
                if not jupyternb:
                    logger.info(f"PubChemPy Error for '{name}': {e}")
                #skip this molecule 
                continue
            #error if general error
            except Exception as e:
                if not jupyternb:
                    logger.info(f"General Error processing '{name}': {e}")
                #skip this molecule
                continue
        
        #add the new category to the additional_categories state variable
        molecules_dict = {item["name"]: item["smiles"] for item in molecules_with_smiles}

        if not jupyternb:
            st.session_state.additional_categories[category_header] = molecules_dict
            
            #store newly created category name temporarily to select it later
            st.session_state.last_created_category = category_header
            
            #immediately update the selected category and set category_is_default to False
            st.session_state.selected_molecule_category = category_header
            st.session_state.category_is_default = False
            
            logger.info(f"Created category: {category_header}, Elements: {molecules_dict}")
            logger.info(f"Updated category: {st.session_state.selected_molecule_category}, isDefault: {st.session_state.category_is_default}")
        
            return "Successfully created category"
        #return category name and molecules
        else:
            return category_header, molecules_dict

    except Exception as e:
        logger.info(f"Error communicating with Gemini or processing response: {e}")
        return f"Error: {e}", []

