import streamlit as st
import time
from sketchem.db.mock_db import get_game, start_game
from sketchem.data.molecules import MOLECULE_CATEGORIES


import logging
logger = logging.getLogger("sketchem_app")

def render_waiting_room():
    """Renders the waiting room for both host and joining players"""

    st.markdown("## Game Lobby")
    st.markdown(f"Player: **{st.session_state.player_name}**")
    st.markdown(f"Game Code: **{st.session_state.game_code}**")

    # Display selected category and molecules if game is created 
    if st.session_state.game_mode == "created_multi":
        st.markdown(f"Selected Category: **{st.session_state.selected_molecule_category}**")
        st.markdown("### Molecules:")
        category = st.session_state.selected_molecule_category[0]
        logger.info(f"Category: {category}")
        
        for molecule in MOLECULE_CATEGORIES[category].keys():
            st.markdown(f"- {molecule}")

    # Get and display current players
    game = get_game(st.session_state.game_code)
    if game:
        st.markdown("### Players:")
        for player_id, player_data in game["players"].items():
            st.markdown(f"- {player_data['name']}")

        

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