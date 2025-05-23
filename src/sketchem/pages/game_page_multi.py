"""
Multiplayer game page for Sketchem.

This file contains the UI and other functions for the multiplayer drawing game.
"""

import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
import time
from streamlit_extras.vertical_slider import vertical_slider
from sketchem.utils.back_button import back_button
from sketchem.db.mock_db import get_game
from sketchem.data.molecules import MOLECULE_CATEGORIES
from streamlit.logger import get_logger
import logging
from sketchem.utils.smiles_validator_ai import get_molecule_with_ai
from sketchem.utils.environment import get_gemini_api_key
from sketchem.db.mock_db import update_player_data
from contextlib import contextmanager

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
</style>
"""

@contextmanager
def st_horizontal():
    """
    Used to reate a horizontal layout for Streamlit elements.
    """
    st.markdown(HORIZONTAL_STYLE, unsafe_allow_html=True)
    with st.container():
        st.markdown('<span class="hide-element horizontal-marker"></span>', unsafe_allow_html=True)
        yield


def save_canvas_as_image(canvas_data):  
    """
    Convert canvas data to image bytes.
    
    Args:
        canvas_data: The canvas data to convert
    Returns:
        The image data as bytes
    """
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None


def toggle_drawing_mode():
    """
    Switch between pen and eraser modes.
    
    Updates the session state to toggle between freedraw and erase modes,
    and restores the last pen color when switching back to pen mode.
    """
    if st.session_state.drawing_mode == "freedraw":
        st.session_state.drawing_mode = "erase"
    else:
        st.session_state.drawing_mode = "freedraw"
        # restore last pen color
        st.session_state.pen_color_selector = st.session_state.last_pen_color

def handle_submission(canvas_result):
    """
    Process a player's drawing submission.
    
    Validates the drawing against the target molecule using AI, updates score
    if correct, and provides feedback to the player.
    
    Args:
        canvas_result: The canvas data containing the player's drawing
    """
    

    
    # Check if canvas is all black (effectively empty)
    img_bytes = save_canvas_as_image(canvas_result.image_data)
    if img_bytes is None or Image.open(io.BytesIO(img_bytes)).getcolors() == [(400*600, (0, 0, 0))]:
        st.session_state.toast_queue = {"message": "Please draw something before submitting!", "icon": "⚠️"}
        st.rerun()
        return
    
    # Get the game to find the target molecule's SMILES
    game = get_game(st.session_state.game_code)
    correct = False
    
    if game and "category" in game:
        category = game["category"]
        target_smiles = None
        
        # Get SMILES from appropriate category
        if category in MOLECULE_CATEGORIES and game.get("category_is_default", True):
            target_smiles = MOLECULE_CATEGORIES[category].get(st.session_state.current_molecule)
        elif not game.get("category_is_default", True) and "additional_categories" in game:
            if category in game["additional_categories"]:
                target_smiles = game["additional_categories"][category].get(st.session_state.current_molecule)
        
        if target_smiles:
            
            logger.info(f"Target smiles: {target_smiles}")

            
            # Validate the drawing against the target SMILES
            api_key = get_gemini_api_key()
            validation_result = get_molecule_with_ai(api_key, img_bytes, target_smiles)

            
            # Handle verification errors
            if isinstance(validation_result, bool):
                correct = validation_result
            elif isinstance(validation_result, str):
                st.session_state.toast_queue = {"message": validation_result, "icon": "❌"}
                st.rerun()
                return
    
    if correct:
        st.session_state.points += 1
        st.session_state.toast_queue = {"message": f"Correct! You drew {st.session_state.current_molecule} correctly.", "icon": "✅"}
        
        # Update score and gameplay_time in database
        try:
            # Calculate elapsed time since game start
            elapsed_time = time.time() - st.session_state.start_time - 5
            
            # Update player data in the database
            
            update_result = update_player_data(elapsed_time)
            
            if not update_result or not update_result.get("success", False):
                logger.error(f"Failed to update player data: {update_result}")
        except Exception as e:
            logger.error(f"Error updating player data: {e}")
        
        # Reset the canvas by incrementing a canvas key counter
        if "canvas_key_counter" not in st.session_state:
            st.session_state.canvas_key_counter = 0
        st.session_state.canvas_key_counter += 1
        
        # Move to next molecule
        select_next_molecule()

        st.rerun()
    else:
        if game.get("hints", False):  # Check if hints are enabled
            if st.session_state.last_gemini_detected_mol:
                try:
                    # Get the first name from the list of names returned by Gemini
                    detected_names = st.session_state.last_gemini_detected_mol.strip().split('\n')
                    if detected_names and detected_names[0].strip():
                        first_name = detected_names[0].strip()
                        st.session_state.toast_queue = {"message": f"Wrong molecule, what you drew looks more like {first_name}", "icon": "☝️"}
                    else:
                        st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
                except Exception as e:
                    logger.error(f"Error processing Gemini response: {e}")
                    st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
            else:
                st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
        else:
            st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
        st.rerun()

def select_next_molecule():
    """
    Select the next molecule for the player to draw.
    
    Gets the next molecule in sequence rather than randomly selecting one.
    If all molecules have been used, marks the player as done.
    """
    # Get the selected category
    game = get_game(st.session_state.game_code)
    if game and "category" in game:  
        category = game["category"]
        if category in MOLECULE_CATEGORIES and game.get("category_is_default", True):
            molecules = list(MOLECULE_CATEGORIES[category].keys())
            if molecules:
                # Get the next molecule in sequence instead of random
                current_index = 0
                if hasattr(st.session_state, 'molecule_index'):
                    current_index = (st.session_state.molecule_index + 1) % len(molecules)
                    # Check if we've gone through all molecules
                    if current_index == 0:
                        st.session_state.player_done = True
                        return
                st.session_state.molecule_index = current_index
                st.session_state.current_molecule = molecules[current_index]
                return
        elif not game.get("category_is_default", True) and "additional_categories" in game:
            if category in game["additional_categories"]:
                molecules = list(game["additional_categories"][category].keys())
                if molecules:
                    # Get the next molecule in sequence instead of random
                    current_index = 0
                    if hasattr(st.session_state, 'molecule_index'):
                        current_index = (st.session_state.molecule_index + 1) % len(molecules)
                        # Check if we've gone through all molecules
                        if current_index == 0:
                            st.session_state.player_done = True
                            return
                    st.session_state.molecule_index = current_index
                    st.session_state.current_molecule = molecules[current_index]
                    return
    return



def handle_skip():
    """
    Handle when a player skips the current molecule.
    
    Resets the canvas, selects the next molecule, and shows a notification.
    """
    # Reset the canvas by incrementing the counter
    if "canvas_key_counter" not in st.session_state:
        st.session_state.canvas_key_counter = 0
    st.session_state.canvas_key_counter += 1
    
    select_next_molecule()
    st.session_state.toast_queue = {"message": "Skipped to next molecule", "icon": "⏭️"}
    st.rerun()

def render_game_page_multi():
    """
    Render the multiplayer game page UI.
    
    Sets up the drawing canvas, color selection buttons, timer, and game state.
    Handles different layouts for mobile and desktop views.
    """
    
    # Try to load the style from file, fall back to inline style if file doesn't exist
    try:
        style_path = '/mount/src/sketchem/src/sketchem/pages/style/multiplayer_game_page_styling.css'
        with open(style_path) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except FileNotFoundError:
        # Fallback inline style
        st.markdown("""
        <style>
        /* White button styling */
        div[data-testid="stButton"] button[key*="White_"] {
            background-color: #ffffff !important;
            border: 2px solid #333 !important;
        }
        
        /* Red button styling */
        div[data-testid="stButton"] button[key*="Red_"] {
            background-color: #ff0000 !important;
            border: 2px solid #333 !important;
        }
        
        /* Blue button styling */
        div[data-testid="stButton"] button[key*="Blue_"] {
            background-color: #0000ff !important;
            border: 2px solid #333 !important;
        }
        
        /* Green button styling */
        div[data-testid="stButton"] button[key*="Green_"] {
            background-color: #00ff00 !important;
            border: 2px solid #333 !important;
        }
        
        /* Yellow button styling */
        div[data-testid="stButton"] button[key*="Yellow_"] {
            background-color: #ffff00 !important;
            border: 2px solid #333 !important;
        }
        
        /* Purple button styling */
        div[data-testid="stButton"] button[key*="Purple_"] {
            background-color: #800080 !important;
            border: 2px solid #333 !important;
        }
        
        /* Selected button styling */
        div[data-testid="stButton"] button[key*="_selected"] {
            border: 3px solid #fff !important; 
            box-shadow: 0 0 0 2px #333 !important;
        }
        
        /* Common button styling */
        div[data-testid="stButton"] button[key*="_color"],
        div[data-testid="stButton"] button[key*="_selected"] {
            width: 28px !important;
            height: 28px !important;
            border-radius: 50% !important;
            padding: 0 !important;
            min-width: unset !important;
            margin: 0 auto !important;
            display: block !important;
        }
        
        /* Eraser button styling */
        div[data-testid="stButton"] button[key="eraser_toggle"] {
            background-color: #f0f0f0 !important;
            border: 2px solid #333 !important;
        }
        </style>
        """, unsafe_allow_html=True)

    st.markdown(
            '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"/>',
            unsafe_allow_html=True,
        )
    
    # initialize more session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"  # default to pen mode
    if "last_pen_color" not in st.session_state:
        st.session_state.last_pen_color = "White"
    if "points" not in st.session_state:
        st.session_state.points = 0
    if "current_molecule" not in st.session_state:
        st.session_state.current_molecule = "Default Molecule"  # Fallback default
        select_next_molecule()
    if "game_over" not in st.session_state:
        st.session_state.game_over = False
    if "start_time" not in st.session_state: #starts timer for player
        st.session_state.start_time = time.time()
    if "player_done" not in st.session_state:
        st.session_state.player_done = False
    
    # Get game info
    game = get_game(st.session_state.game_code)
    game_duration = game.get("game_duration") 
    
    
    padding1, goodcolumn, padding2 = st.columns([1, 3, 1])
    

    with goodcolumn:
        if st.session_state.game_over or st.session_state.player_done:
            back_button(destination=None, label="Leave")
        
        if not (st.session_state.game_over or st.session_state.player_done):

            # Display target molecule
            st.markdown(f"## Please draw: **{st.session_state.current_molecule}**")
            
            # Display game info
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown(f"**Score:** {st.session_state.points}")
            with col2:
                # Get total number of molecules for progress display
                game = get_game(st.session_state.game_code)
                total_molecules = 0
                if game and "category" in game:
                    category = game["category"]
                    if category in MOLECULE_CATEGORIES and game.get("category_is_default", True):
                        total_molecules = len(MOLECULE_CATEGORIES[category])
                    elif not game.get("category_is_default", True) and "additional_categories" in game:
                        if category in game["additional_categories"]:
                            total_molecules = len(game["additional_categories"][category])
                
                # Display progress as current molecule index + 1 out of total
                current_progress = getattr(st.session_state, 'molecule_index', 0) + 1
                st.markdown(f"**Progress:** {current_progress}/{total_molecules}")
            with col3:
                @st.fragment(run_every="1s")
                def timer_fragment(game_duration):
                    elapsed_time = time.time() - st.session_state.start_time
                    remaining_time = max(0, game_duration - elapsed_time)
                    
                    # Check if game is over
                    if remaining_time <= 0 and not st.session_state.game_over:
                        st.session_state.game_over = True
                        st.session_state.toast_queue = {"message": "Game Over!", "icon": "🏁"}
                        st.rerun() #rerun the whole page
                    st.markdown(f"**Time remaining:** {int(remaining_time)}s")
                timer_fragment(game_duration)
            
            # Define color options with hex values
            color_options = {
                "White": "#ffffff",
                "Red": "#ff0000",
                "Blue": "#0000ff",
                "Green": "#00ff00",
                "Yellow": "#ffff00",
                "Purple": "#800080",
            }

            # configure canvas based on mode
            if st.session_state.drawing_mode == "erase":
                current_stroke_color = "#000000"  # Black for eraser on black background
            else:
                # Make sure last_pen_color is a valid key in color_options
                if st.session_state.last_pen_color not in color_options:
                    st.session_state.last_pen_color = "White"  # Default to white if invalid
                current_stroke_color = color_options[st.session_state.last_pen_color]

            
            # Function to handle color selection
            def select_color(color_name):
                st.session_state.last_pen_color = color_name
                st.session_state.pen_color_selector = color_name
                
                # Switch to freedraw mode if currently in eraser mode
                if st.session_state.drawing_mode == "erase":
                    st.session_state.drawing_mode = "freedraw"
            
            # Initialize canvas key counter if it doesn't exist
            if "canvas_key_counter" not in st.session_state:
                st.session_state.canvas_key_counter = 0
            
            # MOBILE LAYOUT
            if st.session_state.is_mobile:
                
                # Update current_stroke_width based on the new size and drawing mode
                if st.session_state.drawing_mode == "erase":
                    current_stroke_width = st.session_state.pen_size + 20
                else:
                    current_stroke_width = st.session_state.pen_size
                
                # Add a centered container for all color buttons
                st.markdown("<div style='display: flex; justify-content: center; width: 100%;'>", unsafe_allow_html=True)
                
                # Create a container for all color buttons
                with st.container():
                    # First row of color buttons
                    with st_horizontal():
                        # White button
                        is_white_selected = st.session_state.last_pen_color == "White"
                        white_key = "White_selected" if is_white_selected else "White_color"
                        st.button("", key=white_key, on_click=select_color, args=("White",), help="Select White")
                        
                        # Red button
                        is_red_selected = st.session_state.last_pen_color == "Red"
                        red_key = "Red_selected" if is_red_selected else "Red_color"
                        st.button("", key=red_key, on_click=select_color, args=("Red",), help="Select Red")
                        
                        # Blue button
                        is_blue_selected = st.session_state.last_pen_color == "Blue"
                        blue_key = "Blue_selected" if is_blue_selected else "Blue_color"
                        st.button("", key=blue_key, on_click=select_color, args=("Blue",), help="Select Blue")
                    
                    # Second row of color buttons
                    with st_horizontal():
                        # Green button
                        is_green_selected = st.session_state.last_pen_color == "Green"
                        green_key = "Green_selected" if is_green_selected else "Green_color"
                        st.button("", key=green_key, on_click=select_color, args=("Green",), help="Select Green")
                        
                        # Yellow button
                        is_yellow_selected = st.session_state.last_pen_color == "Yellow"
                        yellow_key = "Yellow_selected" if is_yellow_selected else "Yellow_color"
                        st.button("", key=yellow_key, on_click=select_color, args=("Yellow",), help="Select Yellow")
                        
                        # Purple button
                        is_purple_selected = st.session_state.last_pen_color == "Purple"
                        purple_key = "Purple_selected" if is_purple_selected else "Purple_color"
                        st.button("", key=purple_key, on_click=select_color, args=("Purple",), help="Select Purple")
                        
                        # Eraser/pen toggle button
                        eraser_key = "pen_toggle" if st.session_state.drawing_mode == "erase" else "eraser_toggle"
                        eraser_help_message = "Switch to pen" if st.session_state.drawing_mode == "erase" else "Switch to eraser"
                        st.button("", on_click=toggle_drawing_mode, key=eraser_key, help=eraser_help_message)
                
                # Close the centered container
                st.markdown("</div>", unsafe_allow_html=True)
                
                try:
                    @st.fragment()
                    def canvas_fragment():
                        canvas_result = st_canvas(
                            stroke_color=current_stroke_color,
                            fill_color="rgba(255, 255, 255, 0)",
                            stroke_width=current_stroke_width if 'current_stroke_width' in locals() else st.session_state.pen_size,
                            background_color="#000000",
                            height=st.session_state.viewport_width,
                            width=st.session_state.viewport_width,  
                            drawing_mode="freedraw",
                            key=f"canvas_{st.session_state.canvas_key_counter}",
                            display_toolbar=True,
                        )
                        return canvas_result
                    canvas_result = canvas_fragment()
                except Exception as e:
                    st.session_state.toast_queue = {"message": f"Canvas error: {e}", "icon": "❌"}
                    st.rerun()
                
                # Horizontal slider below canvas
                st.slider(
                    "Pen Size", 
                    min_value=1, 
                    max_value=20, 
                    value=st.session_state.pen_size,
                    key="pen_size_slider_mobile",
                    on_change=lambda: setattr(st.session_state, "pen_size", st.session_state.pen_size_slider_mobile)
                )
                
            # DESKTOP LAYOUT (original)
            else:
                # Create a centered row of color buttons using columns
                left_spacer, col1, col2, col3, col4, col5, col6, col7, right_spacer = st.columns([1, 1, 1, 1, 1, 1, 1, 1, 1])
                
                with col1:
                    is_white_selected = st.session_state.last_pen_color == "White"
                    white_key = "White_selected" if is_white_selected else "White_color"
                    st.button("", key=white_key, on_click=select_color, args=("White",), help="Select White")
                
                with col2:
                    is_red_selected = st.session_state.last_pen_color == "Red"
                    red_key = "Red_selected" if is_red_selected else "Red_color"
                    st.button("", key=red_key, on_click=select_color, args=("Red",), help="Select Red")
                
                with col3:
                    is_blue_selected = st.session_state.last_pen_color == "Blue"
                    blue_key = "Blue_selected" if is_blue_selected else "Blue_color"
                    st.button("", key=blue_key, on_click=select_color, args=("Blue",), help="Select Blue")
                
                with col4:
                    is_green_selected = st.session_state.last_pen_color == "Green"
                    green_key = "Green_selected" if is_green_selected else "Green_color"
                    st.button("", key=green_key, on_click=select_color, args=("Green",), help="Select Green")
                
                with col5:
                    is_yellow_selected = st.session_state.last_pen_color == "Yellow"
                    yellow_key = "Yellow_selected" if is_yellow_selected else "Yellow_color"
                    st.button("", key=yellow_key, on_click=select_color, args=("Yellow",), help="Select Yellow")
                
                with col6:
                    is_purple_selected = st.session_state.last_pen_color == "Purple"
                    purple_key = "Purple_selected" if is_purple_selected else "Purple_color"
                    st.button("", key=purple_key, on_click=select_color, args=("Purple",), help="Select Purple")
                
                with col7:
                    # Eraser/pen toggle button
                    eraser_key = "pen_toggle" if st.session_state.drawing_mode == "erase" else "eraser_toggle"
                    eraser_help_message = "Switch to pen" if st.session_state.drawing_mode == "erase" else "Switch to eraser"
                    st.button("", on_click=toggle_drawing_mode, key=eraser_key, help=eraser_help_message)
                
                # Canvas and slider layout
                slider_col, canvas_col = st.columns([1, 4])
                
                # Vertical slider on the left
                with slider_col:
                    size = vertical_slider(
                        label="Eraser Size" if st.session_state.drawing_mode == "erase" else "Pen Size",
                        min_value=1,
                        max_value=20,
                        default_value=st.session_state.pen_size,
                        key="pen_size_slider",
                        height=350,
                        track_color="#c0c0c0",  # optional
                        thumb_color="#3a444d",  # optional
                        slider_color="#3a444d", 
                    )
                    
                    # Update pen size in session state
                    if size != st.session_state.pen_size:
                        st.session_state.pen_size = size
                        
                    # Update current_stroke_width based on the new size and drawing mode
                    if st.session_state.drawing_mode == "erase":
                        current_stroke_width = st.session_state.pen_size + 20
                    else:
                        current_stroke_width = st.session_state.pen_size
                
                # Canvas on the right
                with canvas_col:
                    try:
                        @st.fragment()
                        def canvas_fragment():
                            canvas_result = st_canvas(
                                stroke_color=current_stroke_color,
                                fill_color="rgba(255, 255, 255, 0)",
                                stroke_width=current_stroke_width,
                                background_color="#000000",
                                height=400,
                                width=600,
                                drawing_mode="freedraw",
                                key=f"canvas_{st.session_state.canvas_key_counter}",
                                display_toolbar=True,
                            )
                            return canvas_result
                        canvas_result = canvas_fragment()
                    except Exception as e:
                        st.session_state.toast_queue = {"message": f"Canvas error: {e}", "icon": "❌"}
                        st.rerun()
            
            # Buttons row (same for both layouts)
            col1, col2, col3 = st.columns([1, 1, 1])
            
            # Back button
            with col1:
                back_button(destination=None, label="Leave Game")
            
            # Skip button
            with col2:
                if st.button("Skip Molecule", key="skip_btn", use_container_width=True, 
                            disabled=st.session_state.game_over or st.session_state.player_done):
                    handle_skip()
            
            # Submit button
            with col3:
                if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True, 
                            disabled=st.session_state.game_over or st.session_state.player_done):
                    handle_submission(canvas_result)
            
        # Game over screen
        if st.session_state.game_over:
            if game and "category" in game:  
                category = game["category"]
                if category in MOLECULE_CATEGORIES and game.get("category_is_default", True):
                    molecules = list(MOLECULE_CATEGORIES[category].keys())
                elif not game.get("category_is_default", True) and "additional_categories" in game:
                    if category in game["additional_categories"]:
                        molecules = list(game["additional_categories"][category].keys())
            st.markdown(f"## Game Over! Your final score: **{st.session_state.points}/{len(molecules)}**")
            # Hide player_done message when game is over
            st.session_state.player_done = False
        elif st.session_state.player_done:
            if game and "category" in game:  
                category = game["category"]
                if category in MOLECULE_CATEGORIES and game.get("category_is_default", True):
                    molecules = list(MOLECULE_CATEGORIES[category].keys())
                elif not game.get("category_is_default", True) and "additional_categories" in game:
                    if category in game["additional_categories"]:
                        molecules = list(game["additional_categories"][category].keys())
            st.markdown(f"## You're done! Your score: **{st.session_state.points}/{len(molecules)}**")
        
        st.divider()
        
        # Leaderboard
        @st.fragment(run_every="5s")
        def leaderboard_fragment():
            st.markdown("### Leaderboard")
            if game:
                players = game["players"]

                players_list = []
                for player_id, player_data in players.items():
                    players_list.append(player_data)
                # Sort first by score (higher is better), then by time (lower is better)
                sorted_players = sorted(players_list, 
                                    key=lambda x: (x.get("score", 0), -x.get("gameplay_time", 0)), 
                                    reverse=True)
                    
                for i, player in enumerate(sorted_players):
                    player_name = player.get('name', 'Unknown')
                    if player_name == st.session_state.player_name:
                        st.markdown(f"{i+1}. **{player_name} (You)**: {player.get('score', 0)}")
                    else:
                        st.markdown(f"{i+1}. **{player_name}**: {player.get('score', 0)}")
        leaderboard_fragment()           
   

if __name__ == "__main__":
    render_game_page_multi()
