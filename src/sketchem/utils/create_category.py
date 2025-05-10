import streamlit as st
from google import genai
from streamlit.logger import get_logger
import logging
from sketchem.data.molecules import MOLECULE_CATEGORIES
from sketchem.utils.create_category_prompt import create_prompt

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def check_category_is_default(selected_category):
    if selected_category in MOLECULE_CATEGORIES.keys():
        st.session_state.category_is_default = True
    else:
        st.session_state.category_is_default = False

def process_gemini_category_response(response_text):
    """Process Gemini API response and add it to additional_categories"""
    try:
        # Parse the response text into a dictionary
        # Assuming the response is in the format:
        # Category name
        # Molecule1: SMILES1
        # Molecule2: SMILES2
        # ...
        molecules_dict = {}
        lines = response_text.strip().split('\n')
        
        if ":" not in response_text: #Check that gemini's answer is the right formatting; also checks that gemini didn't say it couldn't generate a category as that would likely not contain :
            raise ValueError("Invalid response format: No molecule definitions found (missing ':' separator)")
        
        category_name = "" # Initialize category name so that it doesn't throw an error a few lines down
        for line in lines:
            if ':' not in line:
                category_name = line.strip()
            else:
                molecule, smiles = line.split(':', 1)
                molecules_dict[molecule.strip()] = smiles.strip()
        
        # Add the new category to the additional_categories state var
        st.session_state.additional_categories[category_name] = molecules_dict
        
        # Store the newly created category name temporarily to select it later
        st.session_state.last_created_category = category_name
        
        # Immediately update the selected category and set category_is_default to False
        st.session_state.selected_molecule_category = category_name
        st.session_state.category_is_default = False
        
        logger.info(f"Created category: {category_name}, Elements: {molecules_dict}")
        logger.info(f"Updated category: {st.session_state.selected_molecule_category}, isDefault: {st.session_state.category_is_default}")
        
        return True
    except Exception as e:
        st.error(f"Error processing category: {e}")
        return False

def generate_new_category(api_key, user_prompt):
    """Generate a new molecule category using Gemini AI"""
    # Check for empty API key first
    if not api_key:
        return "Gemini API key not set."
    
    try:
        # Call the Gemini API to get the category
        client = genai.Client(api_key=api_key)
        prompt = create_prompt(user_prompt)
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[prompt],
        )
        
        response_text = response.text.strip()
        
        # Process the response and add to additional_categories
        if process_gemini_category_response(response_text):
            return "Successfully created category"
        else:
            return "Failed to process category"
    
    except Exception as e:
        return f"Gemini API error: {e}"