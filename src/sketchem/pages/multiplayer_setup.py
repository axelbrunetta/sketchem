import streamlit as st
import logging
from sketchem.db.mock_db import create_game, join_game

logger = logging.getLogger("sketchem_app")

# Define molecule categories
MOLECULE_CATEGORIES = {
    "Basic Organic Compounds (1)": {
        "Ethanol": "CCO",
    },
    "Biochemistry (1)": {
        "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
    },
    "Pharmaceuticals (1)": {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    }
}

def handle_join_game(player_name: str, game_code: str):
    """Handles join game button click"""
    if not player_name:
        st.error("Please enter your name first.")
    else:
        with st.spinner("Joining game..."):
            try:
                logger.info(f"Player {player_name} joining game: {game_code}")
                response = join_game(game_code, player_name)
                if response.get("success", False): #Check that joining game worked, defaults to false
                    st.session_state.game_code = game_code
                    st.session_state.player_id = response["player_id"]
                    st.session_state.game_mode = "joined_multi"  #set game mode to join waiting room 
                    st.rerun()
                else:
                    error_msg = response.get("error", "Failed to join game")
                    st.error(error_msg)
            except Exception as e:
                logger.error(f"Error joining game: {e}")
                st.error("Failed to join game")

def handle_create_game(player_name: str):
    """Handles create game button click"""

    with st.spinner("Creating game..."):
        try:
            logger.info(f"Creating new game for player: {player_name}")

            response = create_game(player_name) # This adds the game to the database and returns a dictionary with the game code and the UUID of the player that created it
            #Here we're checking that a game code and player id were actually returned, but also that the whole response returned is =/= to None -> unlikely but could happen if there's a memory issue with the streamlit database 
            if response and "game_code" and "player_id" in response: 
                st.session_state.game_code = response["game_code"]
                st.session_state.player_id = response["player_id"]
                st.session_state.game_mode = "created_multi"
                st.rerun()
            else:
                st.error("Failed to create game")
        except Exception as e:
            logger.error(f"Error creating game: {e}")
            st.error("Failed to create game")


def render_multiplayer_setup():
    """Renders the multiplayer setup page"""

    st.markdown("## Multiplayer Setup")
    
    # Player name input
    player_name = st.text_input("Enter your name:", key="player_name_input")
    
    if not player_name:
        st.warning("Please enter your name to continue")
        return

    
    
    col1, col2 = st.columns(2)
   

    with col1:
        # Use state var show_new_game_button to show/hide button as nothing else worked...
        if "show_new_game_button" not in st.session_state:
            st.session_state.show_create_game_button = True
            
        # Only show the button if the state variable is True
        if st.session_state.show_create_game_button:
            if st.button("New Game", use_container_width=True):
                # Set the state variable to False when clicked
                st.session_state.show_create_game_button = False
                st.rerun()  # Force a rerun to update the UI
        
        if not st.session_state.show_create_game_button:
            # Select game duration
            st.markdown("### Game Settings")
            st.session_state.game_duration = st.slider(
                "Game Duration (seconds)",
                min_value=30,
                max_value=180,
                value=60,
                step=10
            )

            # Molecule category selection
            st.markdown("### Select Molecule Category")
            selected_category = st.selectbox(
                "Choose a category:",
                options=list(MOLECULE_CATEGORIES.keys()),
                key="molecule_category"
            )
            
            if selected_category:
                
                # Display molecules in selected category
                st.markdown(f"**Molecules in {selected_category}:**")
                for molecule in MOLECULE_CATEGORIES[selected_category].keys():
                    st.markdown(f"- {molecule}")

            create_disabled = selected_category is None #Disable button below if no category selected
            if st.button("Create New Game", use_container_width=True, disabled=create_disabled):
                handle_create_game(player_name)

    with col2:
        st.markdown("### Join Existing Game")
        game_code = st.text_input("Enter Game Code:", key="game_code_input").upper()
        
        if st.button("Join Game", use_container_width=True, disabled=not game_code):
            handle_join_game(player_name, game_code)

