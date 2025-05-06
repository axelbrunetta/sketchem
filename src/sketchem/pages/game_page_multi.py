import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
import os
from sketchem.utils.back_button import back_button

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

def render_game_page_multi():
    back_button(destination=None, label="Leave game") #Display back button at the top left
    
    # Load custom CSS for some styling of the page -> gets a little messy with streamlit's light / dark mode stuff
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Need to go up one level to reach the css file
    css_path = os.path.join(current_dir, "styles.css")
    with open(css_path) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    st.markdown("<h1 class='page-title'>Our Page Title Idk</h1>", unsafe_allow_html=True)

    # Initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "canvas_key" not in st.session_state:
        st.session_state.canvas_key = 0

    # Center container for canvas and controls
    with st.container():
        # Create columns for centering
        _, center_col, _ = st.columns([1, 2, 1])
        
        with center_col:
            # Canvas for drawing
            canvas_result = st_canvas(
                fill_color="rgba(255, 255, 255, 0)",  # Transparent fill
                stroke_width=st.session_state.pen_size,
                stroke_color="#ffffff",  # White pen
                background_color="#000000",  # Black background
                height=400,
                width=600,
                drawing_mode="freedraw",
                key=f"canvas_{st.session_state.canvas_key}", # Key used to reset canvas when clear is pressed
                display_toolbar=True,
            )

            # Simple pen size slider
            pen_size = st.slider(
                "Pen Size",
                min_value=1,
                max_value=20,
                value=st.session_state.pen_size, # Need to press canvas twice with pen (i.e. draw two lines) before pen size changes -> need to fix
                key="pen_size_slider"
            )
            st.session_state.pen_size = pen_size

            # Submit button
            if st.button("Submit Drawing", type="primary", key="submit_btn"):
                if canvas_result.image_data is not None:
                    img_bytes = save_canvas_as_image(canvas_result.image_data)
                    if img_bytes:
                        st.success("Drawing submitted successfully!")
                else:
                    st.warning("Please draw something before submitting!")

            # Clear button
            if st.button("Clear Canvas", key="clear_btn"):
                st.session_state.canvas_key += 1
                st.rerun()

if __name__ == "__main__":
    st.set_page_config(
        page_title="Molecular Pictionary",
        layout="centered",
        initial_sidebar_state="collapsed",
    )
    render_game_page()
