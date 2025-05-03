import streamlit as st
from sketchem.db.mock_db import create_game, join_game
from sketchem.data.molecules import MOLECULE_CATEGORIES
from google import genai
from google.genai import types
from streamlit.logger import get_logger
import logging
from sketchem.utils.back_button import back_button

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)






def check_category_is_default(selected_category):
    if selected_category in MOLECULE_CATEGORIES.keys():
        st.session_state.categoryIsDefault = True
    else:
        st.session_state.categoryIsDefault = False


def process_gemini_category_response(response_text):
    """Process Gemini API response and add it to additionalCategories"""
    try:
        # Parse the response text into a dictionary
        # Assuming the response is in the format:
        # Category name
        # Molecule1: SMILES1
        # Molecule2: SMILES2
        # ...
        molecules_dict = {}
        lines = response_text.strip().split('\n')
        
        if ":" not in response_text: #Check that gemini's answer is the right formatting; also checks that gemini didn't say it couldn't generate a category as that would likely not contain :
            raise ValueError("Invalid response format: No molecule definitions found (missing ':' separator)")
        
        category_name = "" # Initialize category name so that it doesn't throw an error a few lines down
        for line in lines:
            if ':' not in line:
                category_name = line.strip()
            else:
                molecule, smiles = line.split(':', 1)
                molecules_dict[molecule.strip()] = smiles.strip()
        
        # Add the new category to the additionalCategories state var
        st.session_state.additionalCategories[category_name] = molecules_dict
        
        # Store the newly created category name temporarily to select it later
        st.session_state.last_created_category = category_name
        
        return True
    except Exception as e:
        st.error(f"Error processing category: {e}")
        return False

def generate_new_category(api_key, user_prompt):
    """Generate a new molecule category using Gemini AI"""
    # Check for empty API key first
    if not api_key:
        return "Gemini API key not set."
    
    try:
        # Call the Gemini API to get the category
        client = genai.Client(api_key=api_key)
        prompt = f"""
Generate a list of molecules that fit most accurately a category described by : "{user_prompt}". 

Please provide 5-10 molecules (except if a number was provided in the "text" from before, in which case use that one for the number of molecules) in the following format:
Category Name (number of molecules)
Molecule 1 Name: SMILES notation
Molecule 2 Name: SMILES notation
...

For example:
Common molecules (3)
Ethanol: CCO
Methane: C
Benzene: C1=CC=CC=C1

⸻

IMPORTANT: Before and fter providing this formatting of the name of the category, name of the molecules and their smiles, do not include ANY other explanations or commentary. Simply output what is asked above.

Be lenient on the category descriptions. If the description if vague try to find molecules related to that description, even if distantly related.

"""
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[prompt],
        )
        
        response_text = response.text.strip()
        
        # Process the response and add to additionalCategories
        if process_gemini_category_response(response_text):
            return "Successfully created category"
        else:
            return "Failed to process category"
    
    except Exception as e:
        return f"Gemini API error: {e}"


def handle_join_game(player_name: str, game_code: str):
    """Handles join game button click"""
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

    back_button(destination=None, label="Back to Home") #Display back button at the top left

    st.markdown("## Multiplayer Setup")
    
    # Player name input
    player_name = st.text_input("Enter your name:", key="player_name_input")
    st.session_state.player_name = player_name 

    if not player_name:
        st.warning("Please enter your name to continue")
        return

    
    
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
            st.button("New Game", on_click=hide_new_game_button, use_container_width=True)
        else:
            # Game settings 
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

            # Function to increment the counter when a new category is added
            def increment_category_counter():
                st.session_state.category_update_counter += 1

            
            #Generate a custom category using gemini
            @st.dialog("Generate a molecule category")
            def openModal():
                st.write(f"What kind of molecule category are you looking for?")
                user_input = st.text_input("")
                if st.button("Submit"):
                    returned_var = generate_new_category(api_key = st.secrets.get("GEMINI_API_KEY", ""), user_prompt = user_input)
                    st.session_state.category_update_counter += 1
                    logger.info(f"Generate category message: {returned_var}")
                    if returned_var == "Successfully created category":
                        st.session_state.toast_queue = {"message": "Successfully created category.", "icon": "✅"}
                    else:
                        st.session_state.toast_queue = {"message": "Failed to create category, try to formulate your query differently.", "icon": "☹️"}
                    st.rerun() #Closes the modal view

            
            
            if st.button("Create a molecule category"):
                openModal()

            

            # Create a new list with all categories
            all_categories = list(MOLECULE_CATEGORIES.keys())
            if hasattr(st.session_state, 'additionalCategories'):
                _ = st.session_state.get("category_update_counter", 0) #forces regeneration of all categories when counter changes
                all_categories.extend(st.session_state.additionalCategories.keys())


            #Select a category
            
            # Use the combined list for the selectbox
            selected_category = st.selectbox( 
                "Choose a category:",
                options=all_categories,
                index=all_categories.index(st.session_state.last_created_category) if hasattr(st.session_state, 'last_created_category') and st.session_state.last_created_category in all_categories else 0,
                key=f"molecule_category_{st.session_state.get('category_update_counter', 0)}"
            ) #Now automatically selects the last created category (when created using AI)
            
            # Clear the last_created_category after using it and update the selected category
            if hasattr(st.session_state, 'last_created_category'):
                # Update the selected category in session state
                st.session_state.selected_molecule_category = st.session_state.last_created_category
                # Check if it's a default category
                check_category_is_default(st.session_state.last_created_category)
                # Clean up
                del st.session_state.last_created_category
                

            if selected_category:
               
                # Display molecules in selected category
                st.session_state.selected_molecule_category = selected_category
                st.markdown(f"**Molecules in {selected_category}:**")
                check_category_is_default(selected_category)
                if st.session_state.categoryIsDefault:
                    for molecule in MOLECULE_CATEGORIES[selected_category].keys(): #display category if default
                        st.markdown(f"- {molecule}")
                else:
                    for molecule in st.session_state.additionalCategories[selected_category].keys(): #display category if ai generated
                        st.markdown(f"- {molecule}")
                

            create_disabled = selected_category is None #Disable button below if no category selected
            if st.button("Create New Game", use_container_width=True, disabled=create_disabled):
                handle_create_game(player_name)

    with col2:
        st.markdown("### Join Existing Game")
        game_code = st.text_input("Enter Game Code:", key="game_code_input").upper()
        
        if st.button("Join Game", use_container_width=True, disabled=not game_code):
            handle_join_game(player_name, game_code)

