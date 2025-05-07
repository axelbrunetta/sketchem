import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
from streamlit_extras.vertical_slider import vertical_slider
from sketchem.utils.environment import is_running_locally
from contextlib import contextmanager

def save_canvas_as_image(canvas_data):  # convert canvas data to png image
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None

# Horizontal layout helper function (copied from waiting_room.py)
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
def st_horizontal(): #Function to create an "inline" block for streamlit elements
    st.markdown(HORIZONTAL_STYLE, unsafe_allow_html=True)
    with st.container():
        st.markdown('<span class="hide-element horizontal-marker"></span>', unsafe_allow_html=True)
        yield

# switch between pen and eraser
def toggle_drawing_mode():
    if st.session_state.drawing_mode == "freedraw":
        st.session_state.drawing_mode = "erase"
    else:
        st.session_state.drawing_mode = "freedraw"
        # restore last pen color
        st.session_state.pen_color_selector = st.session_state.last_pen_color


def render_game_page():
    if is_running_locally():
        import os
        current_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(current_dir, 'style', 'singleplayer_setup_styling.css')) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    else:
        with open('/mount/src/sketchem/src/sketchem/pages/style/singleplayer_setup_styling.css') as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    st.markdown(
            '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"/>',
            unsafe_allow_html=True,
        )
    
    # initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "canvas_key" not in st.session_state:
        st.session_state.canvas_key = 0
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"  # default to pen mode
    if "last_pen_color" not in st.session_state:
        st.session_state.last_pen_color = "White"  # default

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
        current_stroke_color = "#000000"
        current_stroke_width = st.session_state.pen_size + 5
    else:
        current_stroke_color = color_options[st.session_state.last_pen_color]
        current_stroke_width = st.session_state.pen_size

    # Title
    st.markdown("<h1 style='text-align: center; margin-bottom: 20px;'>Single Player Mode</h1>", unsafe_allow_html=True)
    
    # Function to handle color selection
    def select_color(color_name):
        st.session_state.last_pen_color = color_name
        st.session_state.pen_color_selector = color_name
    
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
        if st.session_state.drawing_mode != "erase":
            st.markdown("<div style='height: 50px;'></div>", unsafe_allow_html=True)
            size = vertical_slider(
                label="Pen Size",
                min_value=1,
                max_value=20,
                default_value=st.session_state.pen_size,
                key="pen_size_slider",
                height=300,
            )
            st.session_state.pen_size = size
            current_stroke_width = size
    
    # Canvas on the right
    with canvas_col:
        canvas_result = st_canvas(
            stroke_color=current_stroke_color,
            fill_color="rgba(255, 255, 255, 0)",
            stroke_width=current_stroke_width,
            background_color="#000000",
            height=400,
            width=600,
            drawing_mode=st.session_state.drawing_mode,
            key=f"canvas_{st.session_state.canvas_key}",
            display_toolbar=True,
        )
    
    # Buttons row
    st.markdown("<div style='margin-top: 20px;'></div>", unsafe_allow_html=True)
    submit_col, back_col = st.columns([1, 1])
    
    # Submit button
    with submit_col:
        if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True):
            if canvas_result.image_data is not None:
                img_bytes = save_canvas_as_image(canvas_result.image_data)
                if img_bytes:
                    st.success("Drawing submitted successfully!")
            else:
                st.warning("Please draw something before submitting!")
    
    # Back button
    with back_col:
        if st.button("Back", key="back_btn", use_container_width=True):
            st.session_state.show_back_toast = True
            st.session_state.game_mode = "single_setup"
            st.rerun()


if __name__ == "__main__":
    st.set_page_config(
        page_title="Molecular Pictionary",
        layout="centered",
        initial_sidebar_state="collapsed",
    )
    render_game_page()
