import streamlit as st
from sketchem.utils.environment import is_running_locally

def render_home_page():
    # Center the header
    st.markdown("<h2 style='text-align: center; margin-bottom: 20px;'>Choose Game Mode</h2>", unsafe_allow_html=True)

    # Add custom styling for buttons
    st.markdown("""
    <style>
    div[data-testid="stButton"] > button {
        font-size: 1.1rem;
        font-weight: 500;
    }
    </style>
    """, unsafe_allow_html=True)

    # Create two columns for buttons side by side
    col1, col2 = st.columns(2)

    # Single Player button in left column
    with col1:
        if st.button("Single Player", use_container_width=True):
            st.session_state.game_mode = "single_setup"
            st.rerun()

    # Multiplayer button in right column
    with col2:
        if st.button("Multiplayer", use_container_width=True):
            st.session_state.game_mode = "multiplayer_setup"
            st.rerun()

    # Show info message if running locally (centered)
    if is_running_locally():
        st.info("Note: Multiplayer functionality is limited in local development mode.")