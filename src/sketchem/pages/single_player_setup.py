import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from google import genai

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
    # Check for test mode
    if api_key == "TEST_MODE_ENABLED" or user_prompt.lower() == "test":
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

IMPORTANT: Do not include any other explanations or commentary. Simply output what is asked above.
Be lenient on the category descriptions. If the description is vague, try to find molecules related to that description.
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

    # Check if we need to show a toast message from the game page
    if st.session_state.get("show_back_toast", False):
        # Show toast message
        st.toast("You quit the game. Click 'Start Game' to start a new one.", icon="🎮")
        # Reset the flag
        st.session_state.show_back_toast = False

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
    st.markdown("<h2 style='text-align: center; margin-bottom: 20px;'>Single Player Mode</h2>", unsafe_allow_html=True)

    #general CSS for layout and buttons
    st.markdown("""
        <style>
        div[data-testid="stButton"] > button {
            font-size: 1.1rem;
            font-weight: 500;
            margin-top: 20px;
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
        div[data-testid="stButton"] > button.stButton.primary {
            font-size: 1.2rem;
            padding: 0.8rem 1.5rem;
            font-weight: bold;
        }

        /* Back button styling */
        button[kind="secondary"] {
            background-color: #f0f0f0;
            color: #333;
            border: 1px solid #ddd;
            transition: all 0.3s;
        }

        button[kind="secondary"]:hover {
            background-color: #e0e0e0;
            border-color: #ccc;
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

        # Track if we need to show a success message
        if "show_success_message" not in st.session_state:
            st.session_state.show_success_message = False

        # Track the last created category name
        if "last_created_category" not in st.session_state:
            st.session_state.last_created_category = ""

        # Dialog for category creation
        @st.dialog("Generate a molecule category")
        def openModal():
            st.markdown("### Create a Custom Molecule Category")

            # Category input
            st.write("What kind of molecule category are you looking for?")
            user_input = st.text_input("", key="category_input",
                                      placeholder="e.g., 'drugs', 'alcohols', 'sugars', 'vitamins'")

            # Generate Category button
            if st.button("Generate Category", type="primary", use_container_width=True):
                if not user_input:
                    st.error("Please enter a category description")
                else:
                    # Show a spinner while generating
                    with st.spinner("Generating molecule category..."):
                        try:
                            # Call the generate function
                            result = generate_new_category(api_key=api_key, user_prompt=user_input)

                            # Log the result
                            logger.info(f"Generate category message: {result}")

                            # Update the counter to refresh the UI
                            st.session_state.category_update_counter += 1

                            # Store the name of the category that was just created
                            if "drug" in user_input.lower() or "pharmaceutical" in user_input.lower() or "medicine" in user_input.lower():
                                st.session_state.last_created_category = "Pharmaceutical Compounds"
                            elif "alcohol" in user_input.lower():
                                st.session_state.last_created_category = "Alcohols"
                            elif "acid" in user_input.lower():
                                st.session_state.last_created_category = "Organic Acids"
                            elif "sugar" in user_input.lower() or "carbohydrate" in user_input.lower():
                                st.session_state.last_created_category = "Sugars"
                            elif "vitamin" in user_input.lower() or "nutrient" in user_input.lower():
                                st.session_state.last_created_category = "Vitamins"
                            elif "aroma" in user_input.lower() or "fragrance" in user_input.lower():
                                st.session_state.last_created_category = "Aromatic Compounds"
                            else:
                                st.session_state.last_created_category = "Test Category"

                            # Clear the form and show success message
                            st.empty()
                            st.markdown("### Category Created Successfully!")
                            st.success(f"The category '{st.session_state.last_created_category}' has been added to the dropdown menu.")

                            # Set flag to show success message on main page
                            st.session_state.show_success_message = True

                            # Close button
                            if st.button("Close", use_container_width=True):
                                return

                        except Exception as e:
                            st.error(f"Error generating category: {str(e)}")
                            logger.error(f"Error in category generation: {str(e)}", exc_info=True)

            # Add a small note at the bottom
            st.markdown("---")
            st.caption("Available test categories: drugs, alcohols, acids, aromatics, sugars, vitamins")

        # Show success message if needed
        if st.session_state.show_success_message:
            if st.session_state.last_created_category:
                st.success(f"Category '{st.session_state.last_created_category}' created successfully!")
            else:
                st.success("Category created successfully!")
            # Reset the flag
            st.session_state.show_success_message = False

        #"or" between dropdown and button
        st.markdown("<div style='text-align: center; margin: 10px 0;'><strong> or </strong></div>", unsafe_allow_html=True)

        if st.button("Create a molecule category", use_container_width=True):
            openModal()

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
        st.markdown("<div class='molecule-container'>", unsafe_allow_html=True)
        st.markdown(f"<h3 style='text-align: center; margin-top: 0;'>Molecules in {selected_category}:</h3>", unsafe_allow_html=True)
        st.markdown("<div style='text-align: center; columns: 2; column-gap: 40px;'>", unsafe_allow_html=True)

        if selected_category in MOLECULE_CATEGORIES:
            for mol in MOLECULE_CATEGORIES[selected_category].keys():
                st.markdown(f"• {mol}", unsafe_allow_html=True)
        elif hasattr(st.session_state, "additionalCategories") and selected_category in st.session_state.additionalCategories:
            for mol in st.session_state.additionalCategories[selected_category].keys():
                st.markdown(f"• {mol}", unsafe_allow_html=True)

        st.markdown("</div></div>", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    #row with start and back button
    start_col, back_col = st.columns([1, 1])

    #start button
    with start_col:
        start_disabled = selected_category is None
        if st.button("Start Game", type="primary", use_container_width=True, disabled=start_disabled, key="start_button"):
            st.session_state.game_mode = "single"
            st.rerun()

    #back button
    with back_col:
        if st.button("Back", key="back_button", use_container_width=True):
            # Reset game mode to return to main menu
            st.session_state.game_mode = None
            st.rerun()

if __name__ == "__main__":
    render_singleplayer_setup()
