# The game state is used here to keep track of the parameters like time of game (single or multiplayer), user id, etc -> THIS IS FOR A PER-USER SESSION (i.e. different browser tab, computer, etc)

import streamlit as st

def initialize_game_state():
    """Initialize game state"""

    if "game_mode" not in st.session_state:
        st.session_state.game_mode = None  # None, "single", "multiplayer_setup", "create_multi", "join_multi", "multiplayer", "singleplayer_setup"
    
    if "player_name" not in st.session_state:
        st.session_state.player_name = ""
    
    if "player_id" not in st.session_state:
        st.session_state.player_id = None
        
    if "game_code" not in st.session_state:
        st.session_state.game_code = None
        
    if "game_duration" not in st.session_state:
        st.session_state.game_duration = 60
        
    if "selected_molecule_category" not in st.session_state:
        st.session_state.selected_molecule_category = ""
        
    if "categoryIsDefault" not in st.session_state:
        st.session_state.categoryIsDefault = True

    if "additionalCategories" not in st.session_state:
        st.session_state.additionalCategories = {}

def reset_game_state():
    """Reset game state"""
    pass
    # Will be needed for resetting the various parameters
