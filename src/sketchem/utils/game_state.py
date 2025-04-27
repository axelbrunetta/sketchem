# The game state is used here to keep track of the parameters like time of game (single or multiplayer), user id, etc

import streamlit as st

def initialize_game_state():
    """Initialize game state"""
    if "game_mode" not in st.session_state:
        st.session_state.game_mode = None

def reset_game_state():
    """Reset game state"""
    pass
    # Will be needed for resetting the various parameters


