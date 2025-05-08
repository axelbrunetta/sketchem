import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
import time
import random
from streamlit_extras.vertical_slider import vertical_slider
from sketchem.utils.back_button import back_button
from sketchem.db.mock_db import get_game
from sketchem.data.molecules import MOLECULE_CATEGORIES


@st.fragment()
def timer_fragment(game_duration, run_every="1s"):
    elapsed_time = time.time() - st.session_state.start_time
    remaining_time = max(0, game_duration - elapsed_time)
    
    # Check if game is over
    if remaining_time <= 0 and not st.session_state.game_over:
        st.session_state.game_over = True
        st.session_state.toast_queue = {"message": "Game Over!", "icon": "🏁"}
        st.rerun() #rerun the whole page
    st.markdown(f"**Time remaining:** {int(remaining_time)}s")
    
def save_canvas_as_image(canvas_data):  # convert canvas data to png image
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None


# switch between pen and eraser
def toggle_drawing_mode():
    if st.session_state.drawing_mode == "freedraw":
        st.session_state.drawing_mode = "erase"
    else:
        st.session_state.drawing_mode = "freedraw"
        # restore last pen color
        st.session_state.pen_color_selector = st.session_state.last_pen_color

def handle_submission(canvas_result):
    #MAYBE DO MORE HERE TO CHECK IF IMAGE IS ONLY BLACK??
    if canvas_result.image_data is None:
        st.session_state.toast_queue = {"message": "Please draw something before submitting!", "icon": "⚠️"}
        return
    
  
    correct = True  # need to replace with actual validator
    
    if correct:
        st.session_state.points += 1
        st.session_state.toast_queue = {"message": f"Correct! You drew {st.session_state.current_molecule} correctly.", "icon": "✅"}
        
        #TO DOOO -> Update score in database
        #Also check that score <= nbr of molecules -> if == then win
        
        # Move to next molecule
        select_next_molecule()
        st.rerun()
    else:
        st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}

def select_next_molecule():
    # Get the selected category
    game = get_game(st.session_state.game_code)
    if game and "molecule_category" in game:
        category = game["molecule_category"]
        if category in MOLECULE_CATEGORIES and game.get("category_is_default", True):
            molecules = list(MOLECULE_CATEGORIES[category].keys())
            if molecules:
                # Get the next molecule in sequence instead of random
                current_index = 0
                if hasattr(st.session_state, 'molecule_index'):
                    current_index = st.session_state.molecule_index + 1
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
                        current_index = st.session_state.molecule_index + 1
                    st.session_state.molecule_index = current_index
                    st.session_state.current_molecule = molecules[current_index]
                    return
    



def handle_skip():
    select_next_molecule()
    st.session_state.toast_queue = {"message": "Skipped to next molecule", "icon": "⏭️"}
    st.rerun()

def render_game_page_multi():
    
    with open('/mount/src/sketchem/src/sketchem/pages/style/multiplayer_game_page_styling.css') as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

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
    
    # Get game info
    game = get_game(st.session_state.game_code)
    game_duration = game.get("game_duration") 
    
    
    
    
    # Display game info
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**Score:** {st.session_state.points}")
    with col2:
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

    # Display target molecule
    st.markdown(f"## Please draw: **{st.session_state.current_molecule}** (Category: {game.get('molecule_category', 'None')})")
    
    # Function to handle color selection
    def select_color(color_name):
        st.session_state.last_pen_color = color_name
        st.session_state.pen_color_selector = color_name
        
        # Switch to freedraw mode if currently in eraser mode
        if st.session_state.drawing_mode == "erase":
            st.session_state.drawing_mode = "freedraw"
    
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
        st.markdown("<div style='height: 50px;'></div>", unsafe_allow_html=True)
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
        st.session_state.pen_size = size
        
        # Update current_stroke_width based on the new size and drawing mode
        if st.session_state.drawing_mode == "erase":
            current_stroke_width = st.session_state.pen_size + 20
        else:
            current_stroke_width = st.session_state.pen_size


    # Canvas on the right
    with canvas_col:
        try:
            canvas_result = st_canvas(
                stroke_color=current_stroke_color,
                fill_color="rgba(255, 255, 255, 0)",
                stroke_width=current_stroke_width,
                background_color="#000000",
                height=400,
                width=600,
                drawing_mode="freedraw",
                key=f"canvas",
                display_toolbar=True,
            )
        except Exception as e:
            st.session_state.toast_queue = {"message": f"Canvas error: {e}", "icon": "❌"}
            st.session_state.drawing_mode = "freedraw"
            st.rerun()
    
    # Buttons row
    col1, col2, col3 = st.columns([1, 1, 1])
    
    # Back button
    with col1:
        back_button(destination=None, label="Leave Game")
    
    # Skip button
    with col2:
        if st.button("Skip Molecule", key="skip_btn", use_container_width=True, disabled=st.session_state.game_over):
            handle_skip()
    
    # Submit button
    with col3:
        if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True, disabled=st.session_state.game_over):
            handle_submission(canvas_result)
    
    # Game over screen
    if st.session_state.game_over:
        st.markdown("## Game Over!")
        st.markdown(f"Your final score: **DISPLAY POINTS HERE**")  
        
        # Show leaderboard
        st.markdown("### Leaderboard")
        if game and "players" in game:
            players = game["players"]
            # Sort players by score
            if isinstance(players, list):
                sorted_players = sorted(players, key=lambda x: x.get("score", 0), reverse=True)
                for i, player in enumerate(sorted_players):
                    st.markdown(f"{i+1}. **{player.get('name', 'Unknown')}**: {player.get('score', 0)}") #defaults to 0 unknown if nothing found
                    
   

if __name__ == "__main__":
    render_game_page_multi()
