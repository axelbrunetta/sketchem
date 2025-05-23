
"""
Waiting room page for Sketchem.

This file contains the UI and other functions for the waiting room page.
"""

import streamlit as st
import time
from sketchem.db.mock_db import get_game, start_game
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from contextlib import contextmanager
from sketchem.utils.back_button import back_button
from streamlit_extras.stoggle import stoggle

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

HORIZONTAL_STYLE = """
<style class="hide-element">
    /* Hides the style container and removes the extra spacing */
    .element-container:has(.hide-element) {
        display: none;
    }
    /*
        The selector for >.element-container is necessary to avoid selecting the whole
        body of the streamlit app, which is also a stVerticalBlock.
    */
    div[data-testid="stVerticalBlock"]:has(> .element-container .horizontal-marker) {
        display: flex;
        flex-direction: row !important;
        flex-wrap: wrap;
        gap: 0.5rem;
        align-items: baseline;
    }
    /* Buttons and their parent container all have a width of 704px, which we need to override */
    div[data-testid="stVerticalBlock"]:has(> .element-container .horizontal-marker) div {
        width: max-content !important;
    }
    /* Just an example of how you would style buttons, if desired */
    /*
    div[data-testid="stVerticalBlock"]:has(> .element-container .horizontal-marker) button {
        border-color: red;
    }
    */
</style>
"""
@contextmanager
def st_horizontal(): #Function to create an "inline" block for streamlit elements -> credit: https://gist.github.com/ddorn/decf8f21421728b02b447589e7ec7235
    st.markdown(HORIZONTAL_STYLE, unsafe_allow_html=True)
    with st.container():
        st.markdown('<span class="hide-element horizontal-marker"></span>', unsafe_allow_html=True)
        yield

def render_waiting_room():
    """
    Renders the waiting room for both host and joining players
    
    Displays the game code, player list, category, and molecules. Host can start the game.
    """
    st.empty() #Clears the page -> fix for elements of the multiplayer setup page staying on screen
    padding1, goodcolumn, padding2 = st.columns([1, 3, 1])
    

    with goodcolumn:
        back_button(destination=None, label="Leave game") #Display back button at the top left

        st.markdown("## Game Lobby")
        # Display game code with an option to copy it
        with st_horizontal():
            st.markdown(f"Game Code:") 
            st.code(st.session_state.game_code, language=None)

        game = get_game(st.session_state.game_code)
        if game:
            # Display the original duration (5 seconds less than the actual game duration)
            displayed_duration = game['game_duration'] - 5
            st.markdown(f"Game Duration: **{displayed_duration} seconds**")
            
            # Check if game has started (for both host and joining players)
            if game.get("status") == "active":
                st.session_state.game_mode = "multiplayer"
                st.rerun()
            
        st.divider()

        col3, col4 = st.columns(2)
        with col3:
            st.markdown(f"Your player name: **{st.session_state.player_name}**")

            st.markdown("### Players:")
            for player_id, player_data in game["players"].items():
                if player_data["name"] == st.session_state.player_name:
                    st.markdown(f"- {player_data['name']} (You)")
                else:
                    st.markdown(f"- {player_data['name']}")
            
        with col4:
            if game:
                # Display selected category and molecules for both host and joining players
                category = game["category"]
                st.markdown(f"Selected Category: **{category}**")
                
                # Create molecule list for stoggle
                molecule_list = ""
                if game["category_is_default"]:
                    for molecule in MOLECULE_CATEGORIES[category].keys():
                        molecule_list += f"- {molecule}<br>"
                else:
                    for molecule in game["additional_categories"][category].keys():
                        molecule_list += f"- {molecule}<br>"
                
                # Display molecules using stoggle
                stoggle(
                    f"Molecules in {category}:",
                    f"{molecule_list}",
                )
            
        

        if st.session_state.game_mode == "created_multi": #Only game creator can start the game
            if len(game["players"]) > 1:  
                if st.button("Start Game", type="primary", use_container_width=True): #Button displays if there's at least 2 players
                    start_response = start_game(st.session_state.game_code)
                    if start_response.get("success", False): #Returns the success value for the response, and if there's none defaults to false
                        st.session_state.game_mode = "multiplayer" #Reroute to multiplayer game page
                        st.session_state.start_time = time.time()
                        st.rerun()
                    else:
                        st.error("Failed to start the game")
            else:
                st.info("Waiting for more players to join...")
        else:
            st.info("Waiting for host to start the game...")

        
        time.sleep(2)
        st.rerun()#Needed to auto-refresh the waiting room
