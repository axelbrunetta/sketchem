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

def detect_device_type():
    """Simple device detection based on viewport width"""
    # For now, just return desktop as we'll handle the actual detection in the render function
    return "desktop"

def render_tool_buttons(current_mode):
    container = st.container()
    with container:
        tool_cols = st.columns(2)
        with tool_cols[0]:
            if st.button("✏️", key=f"pen_button_{current_mode}"):
                st.session_state.drawing_mode = "freedraw"
                st.session_state.last_pen_color = "White"
        with tool_cols[1]:
            if st.button("🧽", key=f"eraser_button_{current_mode}"):
                st.session_state.drawing_mode = "freedraw"
                st.session_state.last_pen_color = "Black"

def render_color_buttons(last_pen_color):
    container = st.container()
    with container:
        color_cols = st.columns(len(COLOR_OPTIONS))
        for i, (color_name, color_value) in enumerate(COLOR_OPTIONS.items()):
            with color_cols[i]:
                is_selected = last_pen_color == color_name
                button_key = f"{color_name}_{'selected' if is_selected else 'color'}_{last_pen_color}"
                
                # Add button-specific styling
                st.markdown(f"""
                <style>
                div[data-testid="stButton"] button[key="{button_key}"] {{
                    background-color: {color_value} !important;
                    border: {("3px solid #fff !important; box-shadow: 0 0 0 2px #333;") if is_selected else "2px solid #333 !important"};
                    width: 30px !important;
                    height: 30px !important;
                    border-radius: 50% !important;
                    padding: 0 !important;
                    min-width: unset !important;
                }}
                </style>
                """, unsafe_allow_html=True)
                
                if st.button("", key=button_key, help=f"Select {color_name}"):
                    st.session_state.last_pen_color = color_name
                    st.session_state.drawing_mode = "freedraw"

def render_size_control(last_pen_color, current_pen_size, current_eraser_size):
    container = st.container()
    with container:
        if last_pen_color == "Black":  # Eraser mode
            if st.session_state.device_mode == "desktop":
                eraser_size = vertical_slider(
                    label="Eraser Size",
                    key=f"eraser_size_slider",
                    min_value=5,
                    max_value=20,
                    default_value=current_eraser_size,
                    step=1,
                    thumb_color="#ff4b4b",  # Streamlit red
                    track_color="#0e1117",  # Dark background
                    slider_color="#ff4b4b"  # Streamlit red
                )
                if eraser_size != current_eraser_size:
                    st.session_state.eraser_size = eraser_size
            else:
                eraser_size = st.slider(
                    "Eraser Size",
                    min_value=5,
                    max_value=30,
                    value=current_eraser_size,
                    key=f"eraser_size_slider",
                    on_change=lambda: setattr(st.session_state, 'eraser_size', st.session_state.eraser_size_slider)
                )
            return eraser_size, "eraser"
        else:  # Pen mode
            if st.session_state.device_mode == "desktop":
                pen_size = vertical_slider(
                    label="Pen Size",
                    key=f"pen_size_slider",
                    min_value=1,
                    max_value=20,
                    default_value=current_pen_size,
                    step=1,
                    thumb_color="#ff4b4b",  # Streamlit red
                    track_color="#0e1117",  # Dark background
                    slider_color="#ff4b4b"  # Streamlit red
                )
                if pen_size != current_pen_size:
                    st.session_state.pen_size = pen_size
            else:
                pen_size = st.slider(
                    "Pen Size",
                    min_value=1,
                    max_value=20,
                    value=current_pen_size,
                    key=f"pen_size_slider",
                    on_change=lambda: setattr(st.session_state, 'pen_size', st.session_state.pen_size_slider)
                )
            return pen_size, "pen"

def render_canvas(stroke_color, stroke_width):
    container = st.container()
    with container:
        return st_canvas(
            stroke_color=stroke_color,
            fill_color="rgba(255, 255, 255, 0)",
            stroke_width=stroke_width,
            background_color="#000000",
            height=400,
            width=600,
            drawing_mode="freedraw",
            key="canvas",
            display_toolbar=True,
        )

def save_canvas_as_image(canvas_data):  # convert canvas data to png image
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None

def handle_submission(canvas_result):
    if canvas_result.image_data is None:
        st.session_state.toast_queue = {"message": "Please draw something before submitting!", "icon": "⚠️"}
        return
    
    correct = True  # need to replace with actual validator
    
    if correct:
        st.session_state.points += 1
        st.session_state.toast_queue = {"message": f"Correct! You drew {st.session_state.current_molecule} correctly.", "icon": "✅"}
        select_next_molecule()
        st.rerun()
    else:
        st.session_state.toast_queue = {"message": "Not quite right. Try again!", "icon": "❌"}

def select_next_molecule():
    if "category" not in st.session_state:
        st.error("No category selected. Please go back and choose one.")
        return

    category = st.session_state.category

    if category in MOLECULE_CATEGORIES:
        molecules = list(MOLECULE_CATEGORIES[category].keys())
        if molecules:
            current_index = st.session_state.get("molecule_index", 0)
            if current_index >= len(molecules):
                st.session_state.game_over = True
            else:
                st.session_state.current_molecule = molecules[current_index]
                st.session_state.molecule_index = current_index

def handle_skip():
    select_next_molecule()
    st.session_state.toast_queue = {"message": "Skipped to next molecule", "icon": "⏭️"}
    st.rerun()

def render_game_page():
    # Initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "eraser_size" not in st.session_state:
        st.session_state.eraser_size = 10
    if "canvas_key" not in st.session_state:
        st.session_state.canvas_key = 0
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"
    if "last_pen_color" not in st.session_state:
        st.session_state.last_pen_color = "White"
    if "device_mode" not in st.session_state:
        st.session_state.device_mode = "desktop"
    
    # Add a small toggle for device mode in the corner
    with st.sidebar:
        st.session_state.device_mode = "mobile" if st.checkbox("Mobile Mode", value=st.session_state.device_mode == "mobile") else "desktop"
    
    if "current_molecule" not in st.session_state:
        st.session_state.current_molecule = "Default Molecule"
        select_next_molecule()
    if "game_code" not in st.session_state:
        st.session_state.game_code = "demo_code"
    if "points" not in st.session_state:
        st.session_state.points = 0
    if "molecule_index" not in st.session_state:
        st.session_state.molecule_index = 0
    if "game_over" not in st.session_state:
        st.session_state.game_over = False

    # Layout: controls vs. canvas
    col1, col2 = st.columns([1, 3])

    # Controls column
    with col1:
        # Style for tool buttons
        st.markdown("""
        <style>
        div[data-testid="stButton"] button {
            border-radius: 8px;
            padding: 8px;
            margin: 2px;
            transition: all 0.3s ease;
        }
        div[data-testid="stButton"] button:hover {
            transform: translateY(-2px);
        }
        </style>
        """, unsafe_allow_html=True)

        # Pen and Eraser buttons
        render_tool_buttons(st.session_state.drawing_mode)

        # Highlight active tool
        st.markdown(f"""
        <style>
        div[data-testid="stButton"] button[key="pen_button_{st.session_state.drawing_mode}"] {{
            border: 2px solid {("#4CAF50" if st.session_state.last_pen_color != "Black" else "#ddd")} !important;
        }}
        div[data-testid="stButton"] button[key="eraser_button_{st.session_state.drawing_mode}"] {{
            border: 2px solid {("#ff6b6b" if st.session_state.last_pen_color == "Black" else "#ddd")} !important;
        }}
        </style>
        """, unsafe_allow_html=True)

        # Color buttons (only show when not in eraser mode)
        if st.session_state.last_pen_color != "Black":
            render_color_buttons(st.session_state.last_pen_color)

        st.markdown("<br>", unsafe_allow_html=True)

        # Size slider
        new_size, size_type = render_size_control(
            st.session_state.last_pen_color,
            st.session_state.pen_size,
            st.session_state.eraser_size
        )

    # Canvas column
    with col2:
        # Update stroke color based on current mode
        if st.session_state.last_pen_color == "Black":  # Eraser mode
            current_stroke_color = "#000000"
            current_stroke_width = st.session_state.eraser_size
        else:
            current_stroke_color = COLOR_OPTIONS[st.session_state.last_pen_color]
            current_stroke_width = st.session_state.pen_size

        canvas_result = render_canvas(current_stroke_color, current_stroke_width)

    # Buttons row
    st.markdown("<div style='margin-top: 40px;'></div>", unsafe_allow_html=True)
    back_col, submit_col = st.columns([1, 1])

    with back_col:
        if st.button("Back", key="simple_back_btn", use_container_width=True):
            st.session_state.show_back_toast = True
            st.session_state.game_mode = "single_setup"
            st.rerun()

    with submit_col:
        if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True):
            handle_submission(canvas_result)

    # Game over screen
    if st.session_state.game_over:
        st.markdown("## Game Over!")
        st.markdown(f"Your final score: **{st.session_state.points}**")

if __name__ == "__main__":
    render_game_page()