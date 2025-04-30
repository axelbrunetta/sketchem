import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
import os


def save_canvas_as_image(canvas_data):
    """Convert canvas data to PNG image"""
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])

        # Save to bytes
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None


# Function to toggle drawing mode
def toggle_drawing_mode():
    if st.session_state.drawing_mode == "freedraw":
        st.session_state.drawing_mode = "erase"
    else:
        st.session_state.drawing_mode = "freedraw"

def render_game_page():
    # Load custom CSS for some styling of the page -> gets a little messy with streamlit's light / dark mode stuff
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Need to go up one level to reach the css file
    css_path = os.path.join(current_dir, "styles.css")
    with open(css_path) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    st.markdown("<h1 class='page-title'>Practice individually</h1>", unsafe_allow_html=True)

    # Initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "canvas_key" not in st.session_state:
        st.session_state.canvas_key = 0
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"  # Default to pen mode


    # Center container for canvas and controls
    with st.container():
        # Create columns for centering
        _, center_col, _ = st.columns([1, 2, 1])

        with center_col:

            # Pen size first to debug
            # Let the user choose from a fixed list of colors
            color_options = {
                "White": "#ffffff",
                "Red": "#ff0000",
                "Blue": "#0000ff",
                "Green": "#00ff00",
                "Yellow": "#ffff00",
            }

            # Center the eraser toggle button
            _, center_btn_col, _ = st.columns([1, 2, 1])
            with center_btn_col:
                eraser_label = "Switch to Pen" if st.session_state.drawing_mode == "erase" else "Switch to Eraser"
                st.button(eraser_label, on_click=toggle_drawing_mode, key="eraser_toggle", use_container_width=True)

                # Show a visual indicator of the current mode
                if st.session_state.drawing_mode == "erase":
                    st.markdown("<p style='text-align: center; color: #ff6b6b;'><strong>ERASER MODE</strong></p>", unsafe_allow_html=True)
                else:
                    st.markdown("<p style='text-align: center; color: #4CAF50;'><strong>DRAWING MODE</strong></p>", unsafe_allow_html=True)

            # Only show pen options when not in eraser mode
            if st.session_state.drawing_mode != "erase":
                col1, col2 = st.columns([1, 2])  # Adjust ratio as needed

                with col1:
                    selected_color = st.selectbox("Choose Pen Color", options=list(color_options.keys()), key="pen_color_selector")
                    stroke_color = color_options[selected_color]

                with col2:
                    st.session_state.pen_size = st.slider(
                        "Pen Size",
                        min_value=1,
                        max_value=20,
                        value=st.session_state.pen_size,
                        key="pen_size_slider"
                    )
            else:
                # When in eraser mode, use a default stroke color
                stroke_color = "#ffffff"  # Default color when in eraser mode

        # Configure canvas based on mode
        if st.session_state.drawing_mode == "erase":
            # In erase mode, use background color as stroke color
            current_stroke_color = "#000000"  # Same as background color for erasing
            # Use a larger stroke width for eraser to make it easier to erase
            current_stroke_width = st.session_state.pen_size + 5
        else:
            # In pen mode, use selected color
            current_stroke_color = stroke_color
            current_stroke_width = st.session_state.pen_size

        # Then use that color in the canvas
        # Create a centered container for the canvas
        st.markdown('<div class="canvas-container">', unsafe_allow_html=True)
        canvas_result = st_canvas(
            stroke_color=current_stroke_color,
            fill_color="rgba(255, 255, 255, 0)",  # Transparent fill
            stroke_width=current_stroke_width,
            background_color="#000000",  # Black background
            height=400,
            width=600,
            drawing_mode="freedraw",  # Always use freedraw mode, but change the color for erasing
            key=f"canvas_{st.session_state.canvas_key}",
            display_toolbar=True,
        )
        st.markdown('</div>', unsafe_allow_html=True)

        # Create a row for buttons
        col1, col2 = st.columns(2)

        # Submit button
        with col1:
            if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True):
                if canvas_result.image_data is not None:
                    img_bytes = save_canvas_as_image(canvas_result.image_data)
                    if img_bytes:
                        st.success("Drawing submitted successfully!")
                else:
                    st.warning("Please draw something before submitting!")

        # Clear button
        with col2:
            if st.button("Clear Canvas", key="clear_btn", use_container_width=True):
                st.session_state.canvas_key += 1
                st.rerun()

if __name__ == "__main__":
    st.set_page_config(
        page_title="Molecular Pictionary",
        layout="centered",
        initial_sidebar_state="collapsed",
    )
    render_game_page()
