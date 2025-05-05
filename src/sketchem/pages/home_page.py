import streamlit as st
from sketchem.utils.environment import is_running_locally
import os

def render_home_page():
    """
    Renders the home page with a background image and game mode selection buttons.
    The page has a clean, modern design with a semi-transparent white container
    positioned over a background image.
    """
    # Get the path to the banner image
    current_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    banner_path = os.path.join(current_dir, "banner2.PNG")

    # Convert image to base64 for embedding in CSS
    def get_base64_of_bin_file(bin_file):
        """Convert binary file to base64 string for embedding in CSS"""
        import base64
        with open(bin_file, 'rb') as f:
            return base64.b64encode(f.read()).decode()

    # Main page styling
    st.markdown(f"""
    <style>
    /* Remove default Streamlit padding */
    #root > div:first-child {{ padding-top: 0 !important; }}

    /* Background image styling */
    .stApp {{
        background-image: url("data:image/png;base64,{get_base64_of_bin_file(banner_path)}");
        background-size: contain;
        background-position: center 65%;
        background-repeat: no-repeat;
        background-attachment: fixed;
        min-height: 100vh;
        background-color: #f0f2f6;
        opacity: 0.85;
    }}

    /* White container styling */
    .main .block-container {{
        background-color: rgba(255, 255, 255, 0.85);
        padding: 2rem;
        border-radius: 10px;
        margin: 1rem auto;
        max-width: 800px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        position: relative;
        top: 40vh;
    }}

    /* Button styling */
    div[data-testid="stButton"] > button {{
        font-size: 1.1rem;
        font-weight: 500;
        background-color: rgba(255, 255, 255, 0.9);
        border: 2px solid #FF8C00;
        color: #FF8C00;
        padding: 1rem;
        transition: all 0.3s ease;
    }}

    div[data-testid="stButton"] > button:hover {{
        background-color: #FF8C00;
        color: white;
    }}

    /* ends CSS block */
    </style>
    """, unsafe_allow_html=True)

    # Page content
    st.markdown("<h2 style='text-align: center; margin-bottom: 20px; color: #333;'>Choose Game Mode</h2>",
                unsafe_allow_html=True)

    # Game mode selection buttons
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Single Player", use_container_width=True):
            st.session_state.game_mode = "single_setup"
            st.rerun()
    with col2:
        if not is_running_locally():
            if st.button("Multiplayer", use_container_width=True):
                st.session_state.game_mode = "multiplayer_setup"
                st.rerun()
        else:
            st.info("Multiplayer is only available in the deployed version (Using Streamlit Cloud)")
