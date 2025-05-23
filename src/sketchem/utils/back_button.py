"""
Back button functionality for Sketchem.

This file contains the implementation of the back button that appears
on various pages, handling navigation and cleanup when leaving games.
"""

import streamlit as st
from sketchem.utils.game_state import reset_game_state
from streamlit_extras.stylable_container import stylable_container
from sketchem.db.mock_db import remove_player_from_game


def back_button(destination=None, label="Back", use_container_width=True, key=None):
    """
    Renders a back button that navigates to the specified destination or previous page.
    
    Args:
        destination: The game_mode to return to. If None, will try to determine the previous page.
        label: in case you want to change the label to "Home" or something else
        use_container_width: Whether the button should use the full container width
        key: Optional key for the button
    
    Returns:
        True if the button was clicked, False otherwise
    """
    button_key = key or f"back_button_{destination}" #buttons and elements in streamlit need unique identifier -> prevents errors if multiple back buttons because of a problem with the refresh
    
    #styling
    with stylable_container(
        key="back_button_container",
        css_styles="""
        div[data-testid="stButton"] > button {
            background-color: #88888F !important;
            color: white !important;
            border-radius: 8px !important;
            padding: 0.5em 1em !important;
            width: auto !important;
            border: 2px solid transparent !important;
        }
    """):
        if st.button(label, use_container_width=use_container_width, key=button_key):
            #multiplayer game: remove player from the game
            current_mode = st.session_state.game_mode
            if current_mode in ["created_multi", "joined_multi", "multiplayer"] :
                result = remove_player_from_game(
                    st.session_state.game_code, 
                    st.session_state.player_id
                )
                
                if result.get("success", False):
                    if result.get("game_deleted", False):
                        st.session_state.toast_queue = {"message": "You were the last player. Game was deleted.", "icon": "🗑️"}
                    else:
                        st.session_state.toast_queue = {"message": "Successfully left the game.", "icon": "👋"}
                else:
                    st.session_state.toast_queue = {"message": "Error leaving the game.", "icon": "❌"}
            
            if destination is not None:
                #leads to specific destination
                st.session_state.game_mode = destination
            else:
                #determine previous page based on current game mode
                
                #define reroutes
                navigation_map = { #for now all back buttons go to home page for simplicity
                    "multiplayer": None,  # multiplayer game goes back to home
                    "multiplayer_setup": None,  # setup goes back to home
                    "created_multi": None,  # waiting room (host) goes back home
                    "joined_multi": None,  # waiting room (player) goes back home
                    "guide": None,
                    "singleplayer_setup": None,
                    #need one for single player goes back to home from setup + game page
                }
                
                #previous page set as destination
                st.session_state.game_mode = navigation_map.get(current_mode, None)
            
            #clear any page-specific state if needed
            if current_mode in ["created_multi", "joined_multi", "multiplayer"]:
                reset_game_state()
            st.session_state.selected_molecule_category = None
            #delete session states to set them to default values
            if "pen_size" in st.session_state:
                del st.session_state["pen_size"]
            if "drawing_mode" in st.session_state:
                del st.session_state["drawing_mode"]  
            if "last_pen_color" in st.session_state:
                del st.session_state["last_pen_color"]
            if "points" in st.session_state:
                del st.session_state["points"]
            if "molecule_index" in st.session_state:
                del st.session_state["molecule_index"]
            if "game_over" in st.session_state:
                del st.session_state["game_over"]
            if "start_time" in st.session_state:
                del st.session_state["start_time"]
            if "canvas_key_counter" in st.session_state:
                del st.session_state["canvas_key_counter"]
            if "game_code" in st.session_state:
                del st.session_state["game_code"]
            if "category" in st.session_state:
                del st.session_state["category"]
            if "progress_counter" in st.session_state:
                del st.session_state["progress_counter"]
            if "displayed_molecules" in st.session_state:
                del st.session_state["displayed_molecules"]
            if "current_molecule" in st.session_state:
                del st.session_state["current_molecule"]
            if "molecule_index" in st.session_state:
                del st.session_state["molecule_index"]
                
            st.rerun()
            return True
    
    #return false if button wasn't clicked
    return False