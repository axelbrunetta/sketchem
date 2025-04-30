import streamlit as st
from sketchem.data.molecules import MOLECULE_CATEGORIES
from sketchem.pages.multiplayer_setup import generate_new_category
from streamlit.logger import get_logger
import logging

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def render_singleplayer_setup():
    #get the API key from secrets (if available)
    try:
        api_key = st.secrets["GEMINI_API_KEY"]
    except Exception:
        api_key = None
        logger.warning("Gemini API key not found in secrets")

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
        </style>
    """, unsafe_allow_html=True)

    #columns for setup
    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("### Select Molecule Category")

        all_categories = list(MOLECULE_CATEGORIES.keys())
        update_counter = st.session_state.get("category_update_counter", 0)

        if hasattr(st.session_state, "additionalCategories"):
            all_categories.extend(st.session_state.additionalCategories.keys())

        selected_category = st.selectbox(
            "Choose a category:",
            options=all_categories,
            key=f"molecule_category_{update_counter}"
        )

        @st.dialog("Generate a molecule category")
        def openModal():
            if not api_key:
                st.error("Gemini API key not configured. Custom categories are not available.")
                return
            st.write("What kind of molecule category are you looking for?")
            user_input = st.text_input("")
            if st.button("Submit"):
                result = generate_new_category(api_key=api_key, user_prompt=user_input)
                logger.info(f"Generate category message: {result}")
                st.session_state.category_update_counter += 1
                st.rerun()

        st.markdown("<br>", unsafe_allow_html=True)
        if st.button("Create your own molecule category", use_container_width=True):
            openModal()

    with col2:
        st.markdown("### Game Duration")

        game_duration = st.slider(
            "Time per molecule (seconds):",
            min_value=30,
            max_value=180,
            value=st.session_state.get("game_duration", 60),
            step=10,
            key="single_game_duration"
        )

    #store selections
    if selected_category:
        st.session_state.selected_molecule_category = selected_category
    st.session_state.game_duration = game_duration

    st.markdown("<br>", unsafe_allow_html=True)

    #show molecule list
    if selected_category:
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

    #start button
    _, button_col, _ = st.columns([1, 2, 1])
    with button_col:
        start_disabled = selected_category is None
        if st.button("Start Game", type="primary", use_container_width=True, disabled=start_disabled):
            st.session_state.game_mode = "single"
            st.rerun()