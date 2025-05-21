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
from sketchem.utils.smiles_validator_ai import validate_drawing_with_ai
from sketchem.utils.smiles_validator_ai import get_molecule_with_ai 
from sketchem.utils.environment import get_gemini_api_key
import pubchempy as pcp
from sketchem.db.mock_db import update_player_data
import random

# Initialize logger
logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)
from sketchem.utils.environment import is_running_locally
import os

# Define color options at module level
COLOR_OPTIONS = {
    "White": "#ffffff",
    "Red": "#ff0000",
    "Blue": "#0000ff",
    "Green": "#00ff00",
    "Yellow": "#ffff00",
    "Purple": "#800080",
}

def save_canvas_as_image(canvas_data):  # convert canvas data to png image
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None

def toggle_drawing_mode():
    if st.session_state.drawing_mode == "freedraw":
        st.session_state.drawing_mode = "erase"
    else:
        st.session_state.drawing_mode = "freedraw"
        # restore last pen color
        st.session_state.pen_color_selector = st.session_state.last_pen_color

def handle_submission(canvas_result):
    

    
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
        
        # Add the current molecule to the displayed list
        st.session_state.displayed_molecules.append(st.session_state.current_molecule)

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

        if "progress_counter" not in st.session_state:
            st.session_state.progress_counter = 0
        st.session_state.progress_counter += 1
        
        # Move to next molecule
        select_next_molecule()

        st.rerun()
    else:
        if game.get("hints", False):  # Check if hints are enabled
            if st.session_state.last_gemini_detected_mol:
                try:
                    compounds = pcp.get_compounds(st.session_state.last_gemini_detected_mol, namespace='smiles') 
                    if compounds:
                        compound = compounds[0]
                        st.session_state.toast_queue = {"message": f"Wrong molecule, what you drew looks more like {compound.iupac_name}", "icon": "☝️"}
                    else:
                        st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
                except Exception as e:
                    logger.error(f"PubChem error: {e}")
                    st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
        else:
            st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
        st.rerun()

def select_next_molecule():
    # Get the game to find the category
    game = get_game(st.session_state.game_code)
    if game and "category" in game:
        category = game["category"]
    elif "category" in st.session_state:
        category = st.session_state.category

    if category and (category in MOLECULE_CATEGORIES or category in st.session_state.additional_categories):
        # Get molecules from the appropriate category
        if category in MOLECULE_CATEGORIES:
            molecules = list(MOLECULE_CATEGORIES[category].keys())
        elif category in st.session_state.additional_categories:
            molecules = list(st.session_state.additional_categories[category].keys())
        else:
            molecules = []  # or handle the case where the category is invalid

        if molecules:
            if "displayed_molecules" not in st.session_state:
                st.session_state.displayed_molecules = []

            # Select a random molecule that hasn't been displayed yet
            remaining_molecules = list(set(molecules) - set(st.session_state.displayed_molecules))
            if remaining_molecules:
                st.session_state.current_molecule = random.choice(remaining_molecules)
                st.session_state.displayed_molecules.append(st.session_state.current_molecule)
            else:
                st.session_state.game_over = True  # All molecules have been displayed
        else:
            st.error("Invalid category. Please go back and choose a valid category.")
    else:
        st.error("No category selected. Please go back and choose one.")


def handle_skip():
    # Reset the canvas by incrementing the counter
    if "canvas_key_counter" not in st.session_state:
        st.session_state.canvas_key_counter = 0
    st.session_state.canvas_key_counter += 1
    st.session_state.progress_counter += 1

    select_next_molecule()
    st.session_state.toast_queue = {"message": "Skipped to next molecule", "icon": "⏭️"}
    st.rerun()

def get_molecule_from_category(category):
    """Get molecules from the specified category."""
    # Check if the category is in MOLECULE_CATEGORIES
    if category in MOLECULE_CATEGORIES:
        return list(MOLECULE_CATEGORIES[category].keys())
    elif "additional_categories" in st.session_state and category in st.session_state.additional_categories:
        # If not, check additional categories
        return list(st.session_state.additional_categories[category].keys())
    else:
        return None  # Return None if the category is invalid
    

def render_game_page():
    

    css_path = os.path.join(os.path.dirname(__file__), "style", "singleplayer_game_page_styling.css") if is_running_locally() else '/mount/src/sketchem/src/sketchem/pages/style/singleplayer_game_page_styling.css'
    
    with open(css_path) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    st.markdown(
            '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"/>',
            unsafe_allow_html=True,
        )
       

    # Initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"  # default to pen mode
    if "last_pen_color" not in st.session_state:
        st.session_state.last_pen_color = "White"
    if "points" not in st.session_state:
        st.session_state.points = 0
    if "molecule_index" not in st.session_state:
        st.session_state.molecule_index = 0
    if "game_over" not in st.session_state:
        st.session_state.game_over = False
    if "start_time" not in st.session_state: #starts timer for player
        st.session_state.start_time = time.time()
    if "player_done" not in st.session_state:
        st.session_state.player_done = False
    if "canvas_key_counter" not in st.session_state:
        st.session_state.canvas_key_counter = 0
    if "game_code" not in st.session_state:
        st.session_state.game_code = "demo_code"
    if "category" not in st.session_state:
        st.session_state.category = "Alkanes (8)"  # Replace with an actual category name
    if "progress_counter" not in st.session_state:
        st.session_state.progress_counter = 0

    # Initialize displayed molecules if not already done
    if "displayed_molecules" not in st.session_state:
        st.session_state.displayed_molecules = []

    # Get game info to get category
    game = get_game(st.session_state.game_code)
    if game and "category" in game:
        category = game["category"]
        st.session_state.category = category
        game_duration = game.get("game_duration", 0)  # Set a default value if not found
        
        # Initialize first molecule
        if category is not None:
            molecules = get_molecule_from_category(category)
            if molecules:
                if "current_molecule" not in st.session_state:
                    st.session_state.current_molecule = molecules[0]
                    st.session_state.molecule_index = 0
                    st.session_state.progress_counter = 0  
                    st.session_state.displayed_molecules.append(st.session_state.current_molecule)
    else:
        st.error("No category selected. Please go back and choose one.")
        return
    padding1, goodcolumn, padding2 = st.columns([1, 3, 1])
    

    with goodcolumn:
        
        # Check if the game is over
        if st.session_state.game_over:
            if "category" in st.session_state:
                category = st.session_state.category
                # Get the total number of molecules from the appropriate category
                if category in MOLECULE_CATEGORIES:
                    molecules = list(MOLECULE_CATEGORIES[category].keys())
                elif category in st.session_state.additional_categories:
                    molecules = list(st.session_state.additional_categories[category].keys())
                else:
                    molecules = []  # Handle invalid category case
                st.markdown(f"## Game Over! Your final score: **{st.session_state.points}/{len(molecules)}**")
                
                # Display the back button
                back_button(destination=None, label="Back to Home")
            return  # Exit the function early to prevent rendering other elements

        # Timer fragment logic
        @st.fragment(run_every="1s")
        def timer_fragment():
            elapsed_time = time.time() - st.session_state.start_time
            remaining_time = max(0, game_duration - elapsed_time)
            
            # Check if game is over
            if remaining_time <= 0 and not st.session_state.game_over:
                st.session_state.game_over = True
                st.session_state.toast_queue = {"message": "Game Over!", "icon": "🏁"}
                st.rerun()  # rerun the whole page
            # Remove this line to eliminate the time remaining display at the top left
            # st.markdown(f"**Time remaining:** {int(remaining_time)}s")
        
        timer_fragment()  # Call the timer fragment

        # Create a row of columns for alignment
        padding1, title_column, padding2 = st.columns([1, 2, 1])  # Adjust the proportions as needed

        with title_column:
            # Display target molecule title centered in its column
            st.markdown(f"<h2 style='margin-bottom: 20px; '>Please draw: <strong>{st.session_state.current_molecule}</strong></h2>", unsafe_allow_html=True)

        
        if not (st.session_state.game_over or st.session_state.player_done):
            # Display game info
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown(f"**Score:** {st.session_state.points}")
            with col2:
                st.markdown(f"**Progress:** {st.session_state.progress_counter}/{len(molecules)}")
            with col3:
                @st.fragment(run_every="1s")
                def timer_fragment(game_duration):
                    elapsed_time = time.time() - st.session_state.start_time
                    remaining_time = max(0, game_duration - elapsed_time)
                
                    # Check if game is over
                    if remaining_time <= 0 and not st.session_state.game_over:
                        st.session_state.game_over = True
                        st.session_state.toast_queue = {"message": "Game Over!", "icon": "🏁"}
                        st.rerun()  # rerun the whole page
                    st.markdown(f"**Time remaining:** {int(remaining_time)}s")
                timer_fragment(game_duration)

        # Function to handle color selection
        def select_color(color_name):
            st.session_state.last_pen_color = color_name
            st.session_state.pen_color_selector = color_name
            
            # Switch to freedraw mode if currently in eraser mode
            if st.session_state.drawing_mode == "erase":
                st.session_state.drawing_mode = "freedraw"

        # Create a centered row of color buttons using columns
        css_path = os.path.join(os.path.dirname(__file__), "style", "singleplayer_game_page_styling.css")
        with open(css_path) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

        cols = st.columns(9)
        
        # White button
        with cols[1]:
            is_white_selected = st.session_state.last_pen_color == "White" and st.session_state.drawing_mode != "erase"
            white_key = "White_selected" if is_white_selected else "White_color"
            st.button("", key=white_key, on_click=select_color, args=("White",), help="Select White")
        
        # Red button
        with cols[2]:
            is_red_selected = st.session_state.last_pen_color == "Red" and st.session_state.drawing_mode != "erase"
            red_key = "Red_selected" if is_red_selected else "Red_color"
            st.button("", key=red_key, on_click=select_color, args=("Red",), help="Select Red")
        
        # Blue button
        with cols[3]:
            is_blue_selected = st.session_state.last_pen_color == "Blue" and st.session_state.drawing_mode != "erase"
            blue_key = "Blue_selected" if is_blue_selected else "Blue_color"
            st.button("", key=blue_key, on_click=select_color, args=("Blue",), help="Select Blue")
        
        # Green button
        with cols[4]:
            is_green_selected = st.session_state.last_pen_color == "Green" and st.session_state.drawing_mode != "erase"
            green_key = "Green_selected" if is_green_selected else "Green_color"
            st.button("", key=green_key, on_click=select_color, args=("Green",), help="Select Green")
        
        # Yellow button
        with cols[5]:
            is_yellow_selected = st.session_state.last_pen_color == "Yellow" and st.session_state.drawing_mode != "erase"
            yellow_key = "Yellow_selected" if is_yellow_selected else "Yellow_color"
            st.button("", key=yellow_key, on_click=select_color, args=("Yellow",), help="Select Yellow")
        
        # Purple button
        with cols[6]:
            is_purple_selected = st.session_state.last_pen_color == "Purple" and st.session_state.drawing_mode != "erase"
            purple_key = "Purple_selected" if is_purple_selected else "Purple_color"
            st.button("", key=purple_key, on_click=select_color, args=("Purple",), help="Select Purple")
        
        # Eraser toggle
        with cols[7]:
            eraser_key = "pen_toggle" if st.session_state.drawing_mode == "erase" else "eraser_toggle"
            eraser_help_message = "Switch to pen" if st.session_state.drawing_mode == "erase" else "Switch to eraser"
            st.button("", on_click=toggle_drawing_mode, key=eraser_key, help=eraser_help_message)

        # Configure canvas based on mode
        if st.session_state.drawing_mode == "erase":
            current_stroke_color = "#000000"  # Black for eraser on black background
        else:
            current_stroke_color = COLOR_OPTIONS[st.session_state.last_pen_color]


        # Canvas and slider layout
        slider_col, canvas_col = st.columns([1, 4])
        
        # Vertical slider
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
            # Create a container for the canvas and buttons to ensure same width
            canvas_container = st.container()
            with canvas_container:
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

                # Add spacing before buttons
                st.markdown("<div style='height: 20px'></div>", unsafe_allow_html=True)

        # Create button row outside the nested structure
        button_cols = st.columns([1, 0.1, 1, 0.1, 1])

        # Back button
        with button_cols[0]:
            back_button(destination=None, label="Leave Game")
            
        # Skip button
        with button_cols[2]:
            if st.button("Skip Molecule", key="skip_btn", use_container_width=True, 
                        disabled=st.session_state.game_over):
                handle_skip()

        # Submit button
        with button_cols[4]:
            if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True, 
                        disabled=st.session_state.game_over):
                handle_submission(canvas_result)

if __name__ == "__main__":
    render_game_page()
