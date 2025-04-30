import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from sketchem.pages.multiplayer_setup import generate_new_category
from streamlit.logger import get_logger
import logging

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def render_singleplayer_setup(): #displays setup page for single player mode
    
    # checks if API key is configured
    api_key = None
    try:
        api_key = st.secrets.get("GEMINI_API_KEY")
    except:
        logger.warning("Gemini API key not found in secrets")
    
    # centers the header
    st.markdown("<h2 style='text-align: center; margin-bottom: 20px;'>Single Player Mode</h2>", unsafe_allow_html=True)
    
    # creates a centered column for the content
    _, center_col, _ = st.columns([1, 2, 1])
    
    with center_col:
        # styling buttons
        st.markdown("""
        <style>
        div[data-testid="stButton"] > button {
            font-size: 1.1rem;
            font-weight: 500;
            margin-top: 20px;
        }
        </style>
        """, unsafe_allow_html=True)
        
        st.markdown("### Select Molecule Category")

        #custom created category using gemini
        @st.dialog("Generate a molecule category")
        def openModal():
            if not api_key:
                st.error("Gemini API key not configured. Custom categories are not available.")
                return
                
            st.write(f"What kind of molecule category are you looking for?")
            user_input = st.text_input("")
            if st.button("Submit"):
                returned_var = generate_new_category(api_key=api_key, user_prompt=user_input)
                st.session_state.category_update_counter += 1
                logger.info(f"Generate category message: {returned_var}")
                st.rerun() #Closes the modal view

        
        if st.button("Create your own molecule category"):
            openModal()

        # creates a new list with all categories
        all_categories = list(MOLECULE_CATEGORIES.keys())
        if hasattr(st.session_state, 'additionalCategories'):
            _ = st.session_state.get("category_update_counter", 0) #forces regeneration of all categories when counter changes
            all_categories.extend(st.session_state.additionalCategories.keys())

        # uses the combined list for the selectbox
        selected_category = st.selectbox( 
            "Choose a category:",
            options=all_categories,
            key=f"molecule_category_{st.session_state.get('category_update_counter', 0)}"
        )

        if selected_category:
            # stores the selected category in session state
            st.session_state.selected_molecule_category = selected_category
            
            # shows molecules in selected category 
            st.markdown(f"<h3 style='text-align: center;'>Molecules in {selected_category}:</h3>", unsafe_allow_html=True)
            
            # creates a div for molecule list (div = container used to group/structure content)
            st.markdown("<div style='text-align: center;'>", unsafe_allow_html=True)
            
            # checks whether it's default category or ai-generated
            if selected_category in MOLECULE_CATEGORIES:
                for molecule in MOLECULE_CATEGORIES[selected_category].keys():
                    st.markdown(f"• {molecule}", unsafe_allow_html=True)
            elif hasattr(st.session_state, 'additionalCategories') and selected_category in st.session_state.additionalCategories:
                for molecule in st.session_state.additionalCategories[selected_category].keys():
                    st.markdown(f"• {molecule}", unsafe_allow_html=True)
            
            #closes div
            st.markdown("</div>", unsafe_allow_html=True)

        #game duration slider
        st.markdown("<h3 style='text-align: center;'>Game Duration</h3>", unsafe_allow_html=True)
        
        game_duration = st.slider(
            "Time per molecule (seconds):",
            min_value=30,
            max_value=180,
            value=st.session_state.game_duration,
            step=10,
            key="single_game_duration"
        )
        
        #updates session state with selected duration
        st.session_state.game_duration = game_duration
        
        # start game button
        start_disabled = selected_category is None
        
        if st.button("Start Game", type="primary", use_container_width=True, disabled=start_disabled):
            # switching from setup page to actual game page
            st.session_state.game_mode = "single"
            st.rerun()
