#used to keep track of the parameters like game duration, user id, etc -> THIS IS FOR A PER-USER SESSION (i.e. different browser tab, computer, etc)

import streamlit as st

def initialize_game_state():
    """Initialize game state"""

    if "game_mode" not in st.session_state:
        st.session_state.game_mode = None  #None, "single", "multiplayer_setup", "create_multi", "join_multi", "multiplayer", "singleplayer_setup"

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
        
    if "category_is_default" not in st.session_state:
        st.session_state.category_is_default = True

    if "additional_categories" not in st.session_state:
        st.session_state.additional_categories = {}

    if "returned_category_error" not in st.session_state:
        st.session_state.returned_category_error = False

def reset_game_state():
    """Reset game state"""
    #reset game mode and player information
    st.session_state.game_mode = None
    st.session_state.player_name = ""
    st.session_state.player_id = None
    st.session_state.game_code = None
    
    #reset game settings
    st.session_state.game_duration = 60
    st.session_state.selected_molecule_category = ""
    st.session_state.category_is_default = True
    #keep previously created categories but deselect the current one
    


