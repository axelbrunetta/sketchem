import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
import time
from streamlit_extras.vertical_slider import vertical_slider
from sketchem.utils.back_button import back_button
from sketchem.db.mock_db import get_game
from sketchem.data.molecules import MOLECULE_CATEGORIES

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
    
    correct = True  # need to replace with actual validator
    
    if correct:
        st.session_state.points += 1
        st.session_state.toast_queue = {"message": f"Correct! You drew {st.session_state.current_molecule} correctly.", "icon": "✅"}
        
        # Reset the canvas by incrementing a canvas key counter
        if "canvas_key_counter" not in st.session_state:
            st.session_state.canvas_key_counter = 0
        st.session_state.canvas_key_counter += 1
        
        select_next_molecule()
        st.rerun()
    else:
        st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}
        st.rerun()

def select_next_molecule():
    # Get the game to find the category
    game = get_game(st.session_state.game_code)
    if game and "category" in game:
        category = game["category"]
        if category in MOLECULE_CATEGORIES:
            molecules = list(MOLECULE_CATEGORIES[category].keys())
            if molecules:
                current_index = st.session_state.get("molecule_index", 0)
                next_index = current_index + 1
                if next_index >= len(molecules):
                    st.session_state.game_over = True
                else:
                    st.session_state.current_molecule = molecules[next_index]
                    st.session_state.molecule_index = next_index
                return
            st.error("Invalid category. Please go back and choose a valid category.")
        else:
            st.error("No category selected. Please go back and choose one.")
    else:
        st.error("No category selected. Please go back and choose one.")

def handle_skip():
    # Reset the canvas by incrementing the counter
    if "canvas_key_counter" not in st.session_state:
        st.session_state.canvas_key_counter = 0
    st.session_state.canvas_key_counter += 1
    
    select_next_molecule()
    st.session_state.toast_queue = {"message": "Skipped to next molecule", "icon": "⏭️"}
    st.rerun()

def render_game_page():
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
    if "start_time" not in st.session_state:
        st.session_state.start_time = time.time()
    if "canvas_key_counter" not in st.session_state:
        st.session_state.canvas_key_counter = 0
    if "game_code" not in st.session_state:
        st.session_state.game_code = "demo_code"

    # Get game info to get category
    game = get_game(st.session_state.game_code)
    if game and "category" in game:
        category = game["category"]
        st.session_state.category = category
        
        # Initialize first molecule
        if category in MOLECULE_CATEGORIES:
            molecules = list(MOLECULE_CATEGORIES[category].keys())
            if molecules:
                if "current_molecule" not in st.session_state:
                    st.session_state.current_molecule = molecules[0]
                    st.session_state.molecule_index = 0
    else:
        st.error("No category selected. Please go back and choose one.")
        return

    # Display game info
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**Score:** {st.session_state.points}")
    with col2:
        if "category" in st.session_state:
            molecules = list(MOLECULE_CATEGORIES[st.session_state.category].keys())
            total_molecules = len(molecules)
            st.markdown(f"**Progress:** {st.session_state.molecule_index}/{total_molecules}")

    # Display target molecule
    st.markdown(f"## Please draw: **{st.session_state.current_molecule}**")

    # Function to handle color selection
    def select_color(color_name):
        st.session_state.last_pen_color = color_name
        st.session_state.pen_color_selector = color_name
        
        # Switch to freedraw mode if currently in eraser mode
        if st.session_state.drawing_mode == "erase":
            st.session_state.drawing_mode = "freedraw"

    # Create a centered row of color buttons using columns
    st.markdown("""
        <style>
        [data-testid="column"] {
            width: fit-content !important;
            flex: unset !important;
        }
        [data-testid="column"] > div {
            width: fit-content !important;
        }
        </style>
    """, unsafe_allow_html=True)

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
            
            # Custom CSS to ensure button columns are properly spaced
            st.markdown("""
                <style>
                [data-testid="column"] {
                    padding: 0 !important;
                    margin: 0 !important;
                }
                </style>
            """, unsafe_allow_html=True)
            
            # Create 7 columns for better spacing control
            cols = st.columns([1, 0.1, 1, 0.1, 1])
            
            # Back button
            with cols[0]:
                back_button(destination=None, label="Leave Game")
            
            # Skip button
            with cols[2]:
                if st.button("Skip Molecule", key="skip_btn", use_container_width=True, 
                            disabled=st.session_state.game_over):
                    handle_skip()
            
            # Submit button
            with cols[4]:
                if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True, 
                            disabled=st.session_state.game_over):
                    handle_submission(canvas_result)

    # Game over screen
    if st.session_state.game_over:
        if "category" in st.session_state:
            molecules = list(MOLECULE_CATEGORIES[st.session_state.category].keys())
            st.markdown(f"## Game Over! Your final score: **{st.session_state.points}/{len(molecules)}**")

if __name__ == "__main__":
    render_game_page()