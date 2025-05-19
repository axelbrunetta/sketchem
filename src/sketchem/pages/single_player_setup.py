
import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from sketchem.utils.back_button import back_button
from streamlit_extras.stoggle import stoggle
from sketchem.utils.create_category import get_molecules_for_category_pubchem
import time
from sketchem.utils.environment import is_running_locally, get_gemini_api_key
import os

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def render_singleplayer_setup():
    # Initialize session state variables if they don't exist
    if "category_update_counter" not in st.session_state:
        st.session_state.category_update_counter = 0
    if "additional_categories" not in st.session_state:
        st.session_state.additional_categories = {}
    # Add this line to store the last created category
    if "last_created_category" not in st.session_state:
        st.session_state.last_created_category = None
    if "toast_queue" not in st.session_state:
        st.session_state.toast_queue = None
    # Add this to track if category is default
    if "category_is_default" not in st.session_state:
        st.session_state.category_is_default = True

    # Get the API key using the environment utility function
    api_key = get_gemini_api_key()

    #page title
    st.markdown("<h2 style='margin-bottom: 20px;'>Single Player Setup</h2>", unsafe_allow_html=True)

    css_path = os.path.join(os.path.dirname(__file__), "style", "single_setup_styling.css") if is_running_locally() else '/mount/src/sketchem/src/sketchem/pages/style/single_setup_styling.css'
    
    with open(css_path) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    back_button(destination=None, label="Back to Home")

    #columns for setup
    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("### Molecule Category")

        #prepare categories list but don't display it yet, only when category selected
        all_categories = list(MOLECULE_CATEGORIES.keys())
        update_counter = st.session_state.get("category_update_counter", 0)

        # Add custom categories to the list
        if st.session_state.additional_categories:
            all_categories.extend(st.session_state.additional_categories.keys())

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
        elif selected_category in st.session_state.additional_categories:
            for mol in st.session_state.additional_categories[selected_category].keys():
                molecule_list += f"- {mol}<br>"
        
        stoggle(
            f"Molecules in {selected_category}:",
            f"{molecule_list}",
        )

    st.divider()

    
    if st.button("Start Game", type="secondary", use_container_width=True, key="start_button"):
        if selected_category is None:
            st.toast("Please select a category", icon="⚠️")
        else:
            # Create a game object for single player
            game_code = "single_" + str(int(time.time()))  # Create a unique game code
            st.session_state.game_code = game_code

            # Create game data
            game_data = {
                "code": game_code,
                "status": "active",
                "created_at": int(time.time()),
                "category": selected_category,
                "category_is_default": st.session_state.category_is_default,  # Use the stored value
                "additional_categories": st.session_state.additional_categories,  # Include additional categories
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
