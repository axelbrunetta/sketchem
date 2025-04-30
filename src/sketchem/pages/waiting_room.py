import streamlit as st
import time
from sketchem.db.mock_db import get_game, start_game
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging

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

def st_horizontal(): #Function to create an "inline" block for streamlit elements -> credit: https://gist.github.com/ddorn/decf8f21421728b02b447589e7ec7235
    st.markdown(HORIZONTAL_STYLE, unsafe_allow_html=True)
    with st.container():
        st.markdown('<span class="hide-element horizontal-marker"></span>', unsafe_allow_html=True)
        yield

def render_waiting_room():
    """Renders the waiting room for both host and joining players"""
    st.empty() #Clears the page -> fix for elements of the multiplayer setup page staying on screen

    st.markdown("## Game Lobby")
    st.markdown(f"Your player name: **{st.session_state.player_name}**")

    # Display game code with an option to copy it
    with st_horizontal():
        st.markdown(f"Game Code:") 
        st.code(st.session_state.game_code, language=None)

    st.markdown(f"Game Duration: **{st.session_state.game_duration}**")

    # Display selected category and molecules if game is created 
    if st.session_state.game_mode == "created_multi":
        st.markdown(f"Selected Category: **{st.session_state.selected_molecule_category}**")
        st.markdown("### Molecules:")
        category = st.session_state.selected_molecule_category

        if st.session_state.categoryIsDefault:
            for molecule in MOLECULE_CATEGORIES[category].keys():
                st.markdown(f"- {molecule}")
        else:
            for molecule in st.session_state.additionalCategories[category].keys():
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
