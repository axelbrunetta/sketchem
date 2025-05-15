
import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from google import genai
from streamlit_extras.stoggle import stoggle
from sketchem.utils.create_category import get_molecules_for_category_pubchem
import time

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def process_gemini_category_response(response_text):
    """Process Gemini API response and add it to additionalCategories"""
    try:
        # Parse the response text into a dictionary
        # Assuming the response is in the format:
        # Category name
        # Molecule1: SMILES1
        # Molecule2: SMILES2
        # ...
        molecules_dict = {}
        lines = response_text.strip().split('\n')

        if ":" not in response_text:
            raise ValueError("Invalid response format: No molecule definitions found (missing ':' separator)")

        category_name = None
        for line in lines:
            if ':' not in line:
                category_name = line.strip()
            else:
                molecule, smiles = line.split(':', 1)
                molecules_dict[molecule.strip()] = smiles.strip()

        # Make sure we have a category name
        if not category_name:
            category_name = "Custom Category"

        # Initialize additionalCategories if it doesn't exist
        if "additionalCategories" not in st.session_state:
            st.session_state.additionalCategories = {}

        # Add the new category to the additionalCategories state var
        st.session_state.additionalCategories[category_name] = molecules_dict

        logger.info(f"Added new category: {category_name} with {len(molecules_dict)} molecules")
        return True
    except Exception as e:
        st.error(f"Error processing category: {e}")
        logger.error(f"Error processing category: {e}", exc_info=True)
        return False

def generate_new_category(api_key, user_prompt):
    """Generate a new molecule category using Gemini AI"""
    # Check for test mode - only use test categories when explicitly testing
    if api_key == "TEST_MODE_ENABLED" and (user_prompt.lower() == "test" or 
        any(keyword in user_prompt.lower() for keyword in [
            "alcohol", "acid", "drug", "pharmaceutical", "medicine",
            "aroma", "fragrance", "smell", "sugar", "carbohydrate",
            "vitamin", "nutrient"
        ])):
        logger.info("Using test category")

        # Create a test category based on the prompt
        if "alcohol" in user_prompt.lower():
            test_response = """
Alcohols (5)
Ethanol: CCO
Methanol: CO
Isopropanol: CC(O)C
Butanol: CCCCO
Glycerol: C(C(CO)O)O
"""
        elif "acid" in user_prompt.lower():
            test_response = """
Organic Acids (4)
Acetic Acid: CC(=O)O
Citric Acid: C(C(=O)O)C(CC(=O)O)(C(=O)O)O
Lactic Acid: CC(C(=O)O)O
Formic Acid: C(=O)O
"""
        elif "drug" in user_prompt.lower() or "pharmaceutical" in user_prompt.lower() or "medicine" in user_prompt.lower():
            test_response = """
Pharmaceutical Compounds (8)
Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
Ibuprofen: CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
Paracetamol: CC(=O)NC1=CC=C(C=C1)O
Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Penicillin G: CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
Morphine: CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O
Diazepam: CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C
Fluoxetine: CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F
"""
        elif "aroma" in user_prompt.lower() or "fragrance" in user_prompt.lower() or "smell" in user_prompt.lower():
            test_response = """
Aromatic Compounds (6)
Benzene: C1=CC=CC=C1
Toluene: CC1=CC=CC=C1
Vanillin: COC1=C(C=CC(=C1)C=O)O
Limonene: CC1=CCC(CC1)C(=C)C
Eugenol: COC1=CC(=CC(=C1)O)CC=C
Cinnamaldehyde: C=CC(=O)C=CC1=CC=CC=C1
"""
        elif "sugar" in user_prompt.lower() or "carbohydrate" in user_prompt.lower():
            test_response = """
Sugars (5)
Glucose: C(C1C(C(C(C(O1)O)O)O)O)O
Fructose: C(C(C(C(=O)CO)O)O)O
Sucrose: C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O
Lactose: C(C1C(C(C(C(O1)OC2C(OC(C2O)CO)O)O)O)O)O
Maltose: C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O
"""
        elif "vitamin" in user_prompt.lower() or "nutrient" in user_prompt.lower():
            test_response = """
Vitamins (6)
Vitamin C: C(C(C(C(=O)O)O)O)O
Vitamin A: CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC(=O)C)C)C
Vitamin D3: CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C
Vitamin E: CC1=C(C(=C(C(=C1O)C)C)C)CCC(C)(CCCC(C)CCCC(C)CCCC(C)C)O
Vitamin B1: CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO
Vitamin B12: CC1=C2N3C(=CC4=NC(=C(C5=CC6=NC(=C(C7=CC(=C(N7)C=C8N(C(=C(C9=NC(=C1N2)C(=C9C)C(=O)N)C)C(=O)N)C8CC(=O)N)C)C(=O)N)C6CC(=O)N)C)C(=C4C)C(=O)N)C3CC(=O)N
"""
        else:
            # Default test category
            test_response = """
Test Category (6)
Water: O
Oxygen: O=O
Carbon Dioxide: O=C=O
Methane: C
Ammonia: N
Hydrogen: [H][H]
"""

        # Process the test response
        logger.info(f"Using test response: {test_response}")
        if process_gemini_category_response(test_response.strip()):
            return "Successfully created test category"
        else:
            return "Failed to process test category"

    # Check for empty API key
    if not api_key:
        return "Gemini API key not set."

    try:
        # Call the Gemini API to get the category
        client = genai.Client(api_key=api_key)
        prompt = f"""
Generate a list of molecules that fit most accurately a category described by: "{user_prompt}".

Please provide 5-10 molecules in the following format:
Category Name (number of molecules)
Molecule 1 Name: SMILES notation
Molecule 2 Name: SMILES notation
...

For example:
Common molecules (3)
Ethanol: CCO
Methane: C
Benzene: C1=CC=CC=C1

IMPORTANT: 
1. Do not include any other explanations or commentary. Simply output what is asked above.
2. Be lenient on the category descriptions. If the description is vague, try to find molecules related to that description.
3. For inorganic molecules, include their common names and SMILES notation.
4. For solvents or other chemical categories, provide accurate and relevant molecules.
5. Make sure to include the most common and well-known molecules in the requested category.
"""

        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[prompt],
        )

        response_text = response.text.strip()
        logger.info(f"Gemini API response: {response_text[:100]}...")

        # Process the response and add to additionalCategories
        if process_gemini_category_response(response_text):
            return "Successfully created category"
        else:
            return "Failed to process category"

    except Exception as e:
        logger.error(f"Gemini API error: {e}", exc_info=True)
        return f"Gemini API error: {e}"

def render_singleplayer_setup():
    # Initialize session state variables if they don't exist
    if "category_update_counter" not in st.session_state:
        st.session_state.category_update_counter = 0
    if "additionalCategories" not in st.session_state:
        st.session_state.additionalCategories = {}
    if "toast_queue" not in st.session_state:
        st.session_state.toast_queue = None


    # Get the API key from secrets (if available)
    try:
        api_key = st.secrets["GEMINI_API_KEY"]
        # Check if it's still the placeholder
        if api_key == "your-api-key-here":
            logger.warning("Gemini API key is still the placeholder value")
            api_key = None
    except Exception as e:
        api_key = None
        logger.warning(f"Gemini API key not found in secrets: {e}")

    # Enable test mode for development
    test_mode = True  # Set to True to enable test mode without API key
    if test_mode and api_key is None:
        logger.info("Test mode enabled - using mock API key")
        api_key = "TEST_MODE_ENABLED"

    #page title
    st.markdown("<h2 style='margin-bottom: 20px;'>Single Player Setup</h2>", unsafe_allow_html=True)

    #general CSS for layout and buttons
    st.markdown("""
        <style>
        /* Basic form styling */
        div.stButton > button {
            width: 100%;
            border-radius: 10px;
            font-weight: bold;
            margin-top: 10px;
        }
        
        /* Divider styling */
        hr {
            margin: 20px 0;
            border-color: #f0f0f0;
        }
        
        [data-testid="column"] {
            width: 45% !important;
            padding: 0 2% !important;
        }
        h3 {
            text-align: center !important;
            margin-bottom: 20px !important;
        }
        .stSelectbox, .stSlider {
            margin-top: 20px !important;
        }
        .molecule-container {
            border: 1px solid #ddd;
            border-radius: 10px;
            padding: 15px;
            margin: 20px 0;
            background-color: #f8f9fa;
        }

        /* Create category button styling */
        button[kind="primary"] {
            background: linear-gradient(90deg, #0066cc, #4da6ff, #0066cc) !important;
            background-size: 200% 100% !important;
            color: white !important;
            padding: 15px 32px !important;
            font-size: 16px !important;
            font-weight: bold !important;
            border-radius: 12px !important;
            border: none !important;
            transition: background-position 0.5s ease,
                        transform 0.3s ease, 
                        box-shadow 0.3s ease !important;
        }
        
        button[kind="primary"]:hover {
            background-position: 100% 0 !important;
            transform: scale(1.05) !important;
            box-shadow: 0 0 20px 5px rgba(77, 166, 255, 0.6) !important;
        }

        /* Back and Start Game button styling */
        div[data-testid="stButton"] button[kind="secondary"] {
            background-color: #4f5b66 !important;
            color: white !important;
            padding: 10px 20px !important;
            font-size: 16px !important;
            font-weight: bold !important;
            border-radius: 8px !important;
            border: none !important;
            transition: background-color 0.3s ease !important;
        }

        div[data-testid="stButton"] button[kind="secondary"]:hover {
            background-color: #3a444d !important;
        }
        </style>
    """, unsafe_allow_html=True)

    #columns for setup
    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("### Molecule Category")

        #prepare categories list but don't display it yet, only when category selected
        all_categories = list(MOLECULE_CATEGORIES.keys())
        update_counter = st.session_state.get("category_update_counter", 0)

        if hasattr(st.session_state, "additionalCategories"):
            all_categories.extend(st.session_state.additionalCategories.keys())

        #add a placeholder for the selected category
        if "selected_molecule_category" not in st.session_state:
            st.session_state.selected_molecule_category = None

        #create a selectbox that looks like a dropdown button
        #determine the initial index based on the current selection
        initial_index = 0  #default to "Choose Category"
        if st.session_state.selected_molecule_category in all_categories:
            initial_index = all_categories.index(st.session_state.selected_molecule_category) + 1

        selected_category = st.selectbox(
            "Select a category",  # Add a label for accessibility
            options=["Choose Category"] + all_categories,
            index=initial_index,
            key=f"molecule_category_{update_counter}",
            label_visibility="collapsed"  # Hide the label but keep it for accessibility
        )

        #update the selected category in session state
        if selected_category != "Choose Category":
            st.session_state.selected_molecule_category = selected_category
        else:
            #if "Choose Category" is selected, set to None
            selected_category = None
            st.session_state.selected_molecule_category = None

        # Dialog for category creation
        @st.dialog("Generate a molecule category")
        def openModal():
            st.write(f"What kind of molecule category are you looking for?")
            user_input = st.text_input("", placeholder="e.g., 'drugs', 'alcohols', 'sugars', 'vitamins'")
            if st.button("Submit"):
                # Enable test mode for development
                api_key = "TEST_MODE_ENABLED"
                returned_var = get_molecules_for_category_pubchem(api_key=api_key, user_prompt=user_input)

                st.session_state.category_update_counter += 1

                logger.info(f"Generate category message: {returned_var}")
                
                if returned_var == "Successfully created category":
                    st.session_state.toast_queue = {"message": "Successfully created category.", "icon": "✅"}
                else:
                    st.session_state.toast_queue = {"message": "Failed to create category, try to formulate your query differently.", "icon": "☹️"}
                st.rerun() #Closes the modal view

        #"or" between dropdown and button
        st.markdown("<div style='text-align: center; margin: 10px 0;'><strong> or </strong></div>", unsafe_allow_html=True)

        if st.button("Create a molecule category using AI", key="create_category_button", help="This is an experimental feature, some things may not work as intended.", type="primary", use_container_width=True):
            openModal()

        st.divider()

    with col2:
        st.markdown("### Game Duration (seconds)")

        game_duration = st.slider(
            label="Time per molecule",  #required parameter but will be hidden
            min_value=30,
            max_value=180,
            value=st.session_state.get("game_duration", 60),
            step=10,
            key="single_game_duration",
            label_visibility="collapsed"  #hide label
        )

    #store selections
    if selected_category:
        st.session_state.selected_molecule_category = selected_category
    st.session_state.game_duration = game_duration

    st.markdown("<br>", unsafe_allow_html=True)

    #only show molecule list if a valid category is selected (not "Choose Category")
    if selected_category and selected_category != "Choose Category":
        # Display molecules in selected category
        st.session_state.selected_molecule_category = selected_category
        molecule_list = ""
        if selected_category in MOLECULE_CATEGORIES:
            for mol in MOLECULE_CATEGORIES[selected_category].keys():
                molecule_list += f"- {mol}<br>"
        elif hasattr(st.session_state, "additionalCategories") and selected_category in st.session_state.additionalCategories:
            for mol in st.session_state.additionalCategories[selected_category].keys():
                molecule_list += f"- {mol}<br>"
        
        stoggle(
            f"Molecules in {selected_category}:",
            f"{molecule_list}",
        )

    st.markdown("<br>", unsafe_allow_html=True)

    #row with start and back button
    start_col, back_col = st.columns([1, 1])

    #back button
    with start_col:
        if st.button("Back to Home", key="back_button", use_container_width=True, type="secondary"):
            # Reset game mode to return to main menu
            st.session_state.game_mode = None
            st.rerun()

    #start button
    with back_col:
        start_disabled = selected_category is None
        if st.button("Start Game", type="secondary", use_container_width=True, disabled=start_disabled, key="start_button"):
            # Create a game object for single player
            game_code = "single_" + str(int(time.time()))  # Create a unique game code
            st.session_state.game_code = game_code
            
            # Create game data
            game_data = {
                "code": game_code,
                "status": "active",
                "created_at": int(time.time()),
                "category": selected_category,
                "category_is_default": True,  # Single player always uses default categories for now
                "game_duration": game_duration,
                "hints": False,  # No hints in single player
                "players": {}
            }
            
            # Add game to mock database
            from sketchem.db.mock_db import _games
            _games[game_code] = game_data
            
            st.session_state.game_mode = "single"
            st.rerun()

if __name__ == "__main__":
    render_singleplayer_setup()