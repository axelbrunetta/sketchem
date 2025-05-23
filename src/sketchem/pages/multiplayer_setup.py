
"""
Multiplayer setup page for Sketchem.

This file contains the UI and other functions for the multiplayer setup page.
"""

import streamlit as st
from sketchem.db.mock_db import create_game, join_game, get_game
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from sketchem.utils.back_button import back_button

from sketchem.utils.create_category import check_category_is_default,get_molecules_for_category_pubchem
from streamlit_extras.stoggle import stoggle

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

def handle_join_game(player_name: str, game_code: str):
    """
    Handles join game button click
    
    Args:
        player_name: The name of the player joining the game
        game_code: The code of the game to join
    """
    with st.spinner("Joining game..."):
        try:
            logger.info(f"Player {player_name} joining game: {game_code}")
            
            # First check if the game exists and is not already active
            game = get_game(game_code)
            if game and game.get("status") == "active":
                st.error("This game has already started. Please join another game.")
                return
                
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
    """
    Handles create game button click
    
    Args:
        player_name: The name of the player creating the game
    """

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
    """
    Renders the multiplayer setup page
    
    Displays the player name input, game settings, molecule selection, and buttons for creating and joining games.
    """
    
    
    with open('/mount/src/sketchem/src/sketchem/pages/style/multiplayer_setup_styling.css') as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    padding1, goodcolumn, padding2 = st.columns([1, 2, 1])
    

    with goodcolumn:
        back_button(destination=None, label="Back to Home") #Display back button at the top left
        
        st.markdown("## Multiplayer Setup")
        
        # Player name input
        player_name = st.text_input("Enter your name:", key="player_name_input")
        st.session_state.player_name = player_name 

        if not player_name:
            st.warning("Please enter your name to continue")
            return

        
        st.divider()

        col1, col2 = st.columns(2)
    
    
        if 'category_update_counter' not in st.session_state:
            st.session_state.category_update_counter = 0

        with col1:
            # Initialize state 
            if "show_new_game_button" not in st.session_state:
                st.session_state.show_new_game_button = True
                
            def hide_new_game_button():
                st.session_state.show_new_game_button = False
                
            # Show button or game options based on state 
            if st.session_state.show_new_game_button:
                st.button("New Game", on_click=hide_new_game_button, key="new_game_button", use_container_width=True)
            else:
                # Game settings 
                st.markdown("### Game Settings")

                st.markdown("##### Game duration (seconds)")
                @st.fragment()
                def slider_fragment():
                    
                    displayed_duration = st.slider(
                        "",
                        min_value=180,
                        max_value=600,
                        value=300,
                        step=10
                    )
                    # Add 5 seconds to the actual game duration to account for loading time of the game page -> because of the timer
                    st.session_state.game_duration = displayed_duration + 5
                slider_fragment()

                st.divider()

                st.markdown("##### Hints")
                st.session_state.enable_hints = st.toggle("Enable hints")

                st.divider()

                # Molecule category selection
                st.markdown("##### Molecule Category")


                # Function to increment the counter when a new category is added
                def increment_category_counter():
                    st.session_state.category_update_counter += 1

                
                #Generate a custom category using gemini
                @st.dialog("Generate a molecule category")
                def openModal():
                    st.write(f"What kind of molecule category are you looking for?")
                    user_input = st.text_input("")
                    if st.button("Submit"):
                        
                        returned_var = get_molecules_for_category_pubchem(api_key = st.secrets.get("GEMINI_API_KEY", ""), user_prompt = user_input)

                        st.session_state.category_update_counter += 1

                        logger.info(f"Generate category message: {returned_var}")
                        
                        if returned_var == "Successfully created category":
                            st.session_state.toast_queue = {"message": "Successfully created category.", "icon": "✅"}
                        else:
                            st.session_state.toast_queue = {"message": "Failed to create category, try to formulate your query differently.", "icon": "☹️"}
                        st.rerun() #Closes the modal view

                

                

                # Create a new list with all categories
                all_categories = list(MOLECULE_CATEGORIES.keys())
                if hasattr(st.session_state, 'additional_categories'):
                    _ = st.session_state.get("category_update_counter", 0) #forces regeneration of all categories when counter changes
                    all_categories.extend(st.session_state.additional_categories.keys())


                #Select a category
                
                # Use the combined list for the selectbox
                selected_category = st.selectbox( 
                    "Select a molecule category:",
                    options=all_categories,
                    index=all_categories.index(st.session_state.selected_molecule_category) if  not st.session_state.category_is_default else 0,
                    key=f"molecule_category_{st.session_state.get('category_update_counter', 0)}"
                ) #Now automatically selects the last created category (when created using AI)
                


                
                if st.button("Create a molecule category using AI", key="create_category_button", help="This is an experimental feature, some things may not work as intended."):
                    openModal()

                st.divider()

                
                # Clear the last_created_category after using it and update the selected category
                if hasattr(st.session_state, 'last_created_category'):
                    
                    # Clean up
                    del st.session_state.last_created_category

                

                if selected_category:
                    # Display molecules in selected category
                    st.session_state.selected_molecule_category = selected_category
                    check_category_is_default(selected_category)
                    molecule_list = ""
                    if st.session_state.category_is_default:
                        for molecule in MOLECULE_CATEGORIES[selected_category].keys(): #display category if default
                            molecule_list += f"- {molecule}<br>"
                    else:
                        for molecule in st.session_state.additional_categories[selected_category].keys(): #display category if ai generated
                            molecule_list += f"- {molecule}<br>"
                    
                    stoggle(
                    f"Molecules in {selected_category}:",
                    f"{molecule_list}",
                    )

                st.divider()

                create_disabled = selected_category is None #Disable button below if no category selected
                if st.button("Create New Game", key="create_new_game_button", use_container_width=True, disabled=create_disabled):
                    handle_create_game(player_name)

    with col2:
        st.markdown("### Join Existing Game")
        game_code = st.text_input("Enter Game Code:", key="game_code_input").upper()
        
        if st.button("Join Game", key="join_game_button", use_container_width=True, disabled=not game_code):
            handle_join_game(player_name, game_code)

