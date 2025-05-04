import streamlit as st
from sketchem.db.mock_db import create_game, join_game
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from sketchem.utils.back_button import back_button
from sketchem.utils.create_category import check_category_is_default, generate_new_category
import os
from streamlit_extras.stylable_container import stylable_container

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

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

    
    with open('/mount/src/sketchem/src/sketchem/pages/style/multiplayer_setup_styling.css') as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
        
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
            
            st.session_state.enable_hints = st.toggle("Enable hints")

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

            
            with stylable_container(
            key="join_game_container",
            css_styles="""
                *,
*:after,
*:before {
	box-sizing: border-box;
}

:root {
	--transition: 0.25s;
	--spark: 1.8s;
	/* --active is now controlled purely by :has() or button state */
}

body {
	--active: 0; /* Default state */
	background: hsl(
		260
		calc(var(--active) * 97%) /* Only becomes purple when --active is 1 */
		6%
	);
	display: grid;
	place-items: center;
	min-height: 100vh;
	font-family: system-ui, sans-serif;
	transition: background var(--transition);
	overflow: hidden; /* Keep to prevent scrollbars if button scales near edge */
}

/* --- Use :has() to activate body background --- */
/* Supported in modern browsers */
@supports(selector(:has(:is(+ *)))) {
	body:has(.sparkle-button button:is(:hover, :focus-visible)) {
		--active: 1; /* Activate background change */
	}
}
/* Fallback for browsers without :has() - body won't change color */


button {
	--cut: 0.1em;
	/* --active will be 1 on hover/focus */
	--active: 0;
	--bg:
		radial-gradient(
			40% 50% at center 100%,
			hsl(270 calc(var(--active) * 97%) 72% / var(--active)), /* Only visible when --active is 1 */
			transparent
		),
		radial-gradient(
			80% 100% at center 120%,
			hsl(260 calc(var(--active) * 97%) 70% / var(--active)), /* Only visible when --active is 1 */
			transparent
		),
		hsl(260 calc(var(--active) * 97%) calc((var(--active) * 44%) + 12%)); /* Base color changes */
	background: var(--bg);
	font-size: 2rem;
	font-weight: 500;
	color: hsl(0 0% 95%); /* Set a base text color */
	border: 0;
	cursor: pointer;
	padding: 0.9em 1.3em;
	display: flex;
	align-items: center;
	gap: 0.25em;
	white-space: nowrap;
	border-radius: 100px;
	position: relative;
	box-shadow:
		0 0 calc(var(--active) * 6em) calc(var(--active) * 3em) hsl(260 97% 61% / 0.75), /* Glow on active */
		0 0.05em 0 0 hsl(260 calc(var(--active) * 97%) calc((var(--active) * 50%) + 30%)) inset, /* Inner highlight */
		0 -0.05em 0 0 hsl(260 calc(var(--active) * 97%) calc(var(--active) * 60%)) inset; /* Inner shadow */
	transition: box-shadow var(--transition), scale var(--transition), background var(--transition);
	scale: calc(1 + (var(--active) * 0.1)); /* Scale up on active */
}

/* Set --active directly on the button for its internal styles */
button:is(:hover, :focus-visible) {
	--active: 1;
}

button:active {
  /* Slightly reduce scale on click for feedback */
  scale: calc(1 + (var(--active) * 0.08));
}

svg {
	overflow: visible !important; /* Important for bounce effect */
}


/* Trigger bounce only when button is hovered/focused */
button:is(:hover, :focus-visible) .sparkle path {
	animation-name: bounce;
}

@keyframes bounce {
	35%, 65% {
		scale: var(--scale);
	}
}



/* --- Text Styling --- */
.text {
	translate: 2% -6%; /* Slight adjustment */
	letter-spacing: 0.01ch;
    /* Gradient changes based on button's --active state */
	background: linear-gradient(90deg,
        hsl(0 0% calc((var(--active) * 100%) + 65%)), /* Brightens on active */
        hsl(0 0% calc((var(--active) * 100%) + 26%))  /* Brightens on active */
    );
	-webkit-background-clip: text;
	background-clip: text;
	color: transparent;
	transition: background var(--transition);
    z-index: 2; /* Ensure text is above backdrop */
    position: relative; /* Needed for z-index */
}


/* --- Main Container (Optional but good practice) --- */
.sparkle-button {
	position: relative; /* Establishes stacking context */
}
            """
            ):
                if st.button("Create a molecule category using AI"):
                    openModal()

            

            # Create a new list with all categories
            all_categories = list(MOLECULE_CATEGORIES.keys())
            if hasattr(st.session_state, 'additional_categories'):
                _ = st.session_state.get("category_update_counter", 0) #forces regeneration of all categories when counter changes
                all_categories.extend(st.session_state.additional_categories.keys())


            #Select a category
            
            # Use the combined list for the selectbox
            selected_category = st.selectbox( 
                "Choose a category:",
                options=all_categories,
                index=all_categories.index(st.session_state.selected_molecule_category) if  not st.session_state.category_is_default else 0,
                key=f"molecule_category_{st.session_state.get('category_update_counter', 0)}"
            ) #Now automatically selects the last created category (when created using AI)
            
            # Clear the last_created_category after using it and update the selected category
            if hasattr(st.session_state, 'last_created_category'):
                
                # Clean up
                del st.session_state.last_created_category
                

            if selected_category:
               
                # Display molecules in selected category
                st.session_state.selected_molecule_category = selected_category
                st.markdown(f"**Molecules in {selected_category}:**")
                check_category_is_default(selected_category)
                if st.session_state.category_is_default:
                    for molecule in MOLECULE_CATEGORIES[selected_category].keys(): #display category if default
                        st.markdown(f"- {molecule}")
                else:
                    for molecule in st.session_state.additional_categories[selected_category].keys(): #display category if ai generated
                        st.markdown(f"- {molecule}")
                

            create_disabled = selected_category is None #Disable button below if no category selected
            if st.button("Create New Game", use_container_width=True, disabled=create_disabled):
                handle_create_game(player_name)

    with col2:
        st.markdown("### Join Existing Game")
        game_code = st.text_input("Enter Game Code:", key="game_code_input").upper()
        
        if st.button("Join Game", use_container_width=True, disabled=not game_code):
            handle_join_game(player_name, game_code)

