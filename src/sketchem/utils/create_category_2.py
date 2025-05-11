from google import genai
import pubchempy as pcp
from streamlit.logger import get_logger
import logging
import streamlit as st

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)


def get_molecules_for_category_pubchem(api_key, user_prompt):
    full_prompt = f"""
Generate a list of molecule names that fit most accurately a category described by : "{user_prompt}".

Please provide 5-10 molecule names (except if a number was provided in the "text" from before, in which case use that one for the number of molecules).

Output Format:
Provide your response as a simple list of molecule names, each on a new line. Do NOT include SMILES notation, category names, or any other text.

Category Name
Molecule 1 Name
Molecule 2 Name
...

For example, if the category was "Common Alcohols":
Common Alcohols
Methanol
Ethanol
Isopropanol

Or if the category was "Halogenated Solvents":
Halogenated Solvents
Chloroform
Dichloromethane
⸻

IMPORTANT Rules to Follow for Molecule Name Selection and Output:

1.  Names and Category Only: Your entire output should be the category name (with number of molecules) followed by a list of molecule names, each on a new line. Do not include ANY other explanations, commentary, or SMILES strings.
2.  Real Molecules Only: DO NOT invent molecules or use molecules that do not exist in real life (e.g., do NOT use ones that exist in movies or stories, or are fictional / imaginary). Use common, well-established chemical names.
3.  Simple Structures: Focus on molecules that can be easily represented by simple chemical structures (suitable for hand drawing). Avoid names that refer to:
    *   DNA, RNA
    *   Graphite, Diamond (extended covalent networks)
    *   Alloys
    *   Fullerenes, Carbon nanotubes
    *   Complex / Very Large Proteins or Polymers
4.  Molecules, Not Atoms: Do NOT output names of elemental atoms only. Do NOT output things like:
    *   Osmium
    *   Bismuth
    *   Radon
    *   Krypton
    *   Iron
    *   Gold
5.  Category Interpretation: Be lenient on the category descriptions. If the description is vague, try to find molecule names related to that description (even if distantly related) while still respecting all aforementioned rules.
6.  Clarity of Names: Try to use names that are likely to be found in chemical databases like PubChem (e.g., common names or IUPAC names). Avoid overly obscure or ambiguous names.
"""

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

        logger.info(f"Created category: {category_header}, Elements: {molecule_names_from_ai}")

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
                        logger.info(f"Found SMILES for '{name}': {compound.canonical_smiles}")
                    else:
                        logger.info(f"Warning: No canonical SMILES found for '{name}' in PubChem.")
                        
                        return "NO_SMILES_FOUND_IN_PUBCHEM"
                else:
                    logger.info(f"Warning: No molecule found for '{name}' in PubChem.")
                    return "NAME_NOT_FOUND_IN_PUBCHEM"
            except pcp.PubChemPyError as e:
                logger.info(f"PubChemPy Error for '{name}': {e}")
                return "PUBCHEM_API_ERROR"
            except Exception as e:
                logger.info(f"General Error processing '{name}': {e}")
                return "PROCESSING_ERROR"
        
        # Add the new category to the additional_categories state var
        molecules_dict = {item["name"]: item["smiles"] for item in molecules_with_smiles}

        st.session_state.additional_categories[category_header] = molecules_dict
        
        # Store the newly created category name temporarily to select it later
        st.session_state.last_created_category = category_header
        
        # Immediately update the selected category and set category_is_default to False
        st.session_state.selected_molecule_category = category_header
        st.session_state.category_is_default = False
        
        logger.info(f"Created category: {category_header}, Elements: {molecules_dict}")
        logger.info(f"Updated category: {st.session_state.selected_molecule_category}, isDefault: {st.session_state.category_is_default}")
        
        return "Successfully created category"

    except Exception as e:
        logger.info(f"Error communicating with Gemini or processing response: {e}")
        return f"Error: {e}", []

