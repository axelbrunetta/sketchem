import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io


def save_canvas_as_image(canvas_data): #convert canvas data to png image
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None


#switch between pen and eraser
def toggle_drawing_mode():
    if st.session_state.drawing_mode == "freedraw":
        st.session_state.drawing_mode = "erase"
    else:
        st.session_state.drawing_mode = "freedraw"
        #restore last pen color
        st.session_state.pen_color_selector = st.session_state.last_pen_color



def render_game_page():
    #initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "canvas_key" not in st.session_state:
        st.session_state.canvas_key = 0
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"  #default to pen mode
    if "last_pen_color" not in st.session_state:
        st.session_state.last_pen_color = "#ffffff"  #default


    #color options
    color_options = {
        "White": "#ffffff",
        "Red": "#ff0000",
        "Blue": "#0000ff",
        "Green": "#00ff00",
        "Yellow": "#ffff00",
    }

    #configure canvas based on mode
    if st.session_state.drawing_mode == "erase":
        current_stroke_color = "#000000"
        current_stroke_width = st.session_state.pen_size + 5
        stroke_color = "#ffffff"
    else:
        stroke_color = color_options["White"]
        current_stroke_color = stroke_color
        current_stroke_width = st.session_state.pen_size

    #css for spacing and layout
    st.markdown("""
    <style>
    /* Add top margin to push content lower */
    .main > .block-container {
        padding-top: 80px !important;
    }

    h1 {
        text-align: center;
        color: #4CAF50;
        margin-bottom: 0 !important;
    }

    [data-testid="stVerticalBlock"] {
        gap: 0 !important;
    }

    div[data-testid="stButton"] > button {
        font-size: 1.2rem;
        padding: 0.8rem 1.5rem;
        font-weight: bold;
    }
    </style>
    """, unsafe_allow_html=True)

    #title
    st.title("Single Player Mode")
    st.markdown("""
    <div style='margin-top: 80px;'></div>
    """, unsafe_allow_html=True)


    #layout: controls vs. canvas
    col1, col2 = st.columns([1, 3])

    #controls
    with col1:
        eraser_label = "Switch to Pen" if st.session_state.drawing_mode == "erase" else "Switch to Eraser"
        st.button(eraser_label, on_click=toggle_drawing_mode, key="eraser_toggle", use_container_width=True)

        if st.session_state.drawing_mode == "erase":
            st.markdown("<p style='text-align: center; color: #ff6b6b;'><strong>ERASER MODE</strong></p>", unsafe_allow_html=True)
        else:
            st.markdown("<p style='text-align: center; color: #4CAF50;'><strong>DRAWING MODE</strong></p>", unsafe_allow_html=True)

        st.markdown("<br><br><br>", unsafe_allow_html=True)

        if st.session_state.drawing_mode != "erase":
            selected_color = st.selectbox("Choose Pen Color", options=list(color_options.keys()), key="pen_color_selector")
            stroke_color = color_options[selected_color]
            st.session_state.last_pen_color = selected_color  #stores last pen color (no default color after having used eraser)
            current_stroke_color = stroke_color

        st.markdown("<br><br><br><br><br>", unsafe_allow_html=True)

        if st.session_state.drawing_mode != "erase":
            st.session_state.pen_size = st.slider(
                "Pen Size",
                min_value=1,
                max_value=20,
                value=st.session_state.pen_size,
                key="pen_size_slider"
            )
            current_stroke_width = st.session_state.pen_size

    #canvas
    with col2:
        canvas_result = st_canvas(
            stroke_color=current_stroke_color,
            fill_color="rgba(255, 255, 255, 0)",
            stroke_width=current_stroke_width,
            background_color="#000000",
            height=400,
            width=600,
            drawing_mode="freedraw",
            key=f"canvas_{st.session_state.canvas_key}",
            display_toolbar=True,
        )

    #buttons row
    st.markdown("<br>", unsafe_allow_html=True)

    #add custom styling for buttons
    st.markdown("""
    <style>
    /* Make buttons more prominent */
    div[data-testid="stButton"] > button {
        font-size: 1.1rem;
        font-weight: 500;
    }
    </style>
    """, unsafe_allow_html=True)

    #2 equal columns for buttons
    submit_col, back_col = st.columns([1, 1])

    #submit button (left)
    with submit_col:
        if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True):
            if canvas_result.image_data is not None:
                img_bytes = save_canvas_as_image(canvas_result.image_data)
                if img_bytes:
                    st.success("Drawing submitted successfully!")
            else:
                st.warning("Please draw something before submitting!")

    #back button (right)
    with back_col:
        if st.button("Back", key="back_btn", use_container_width=True):
            #store a flag in session state to show the toast on the next page
            st.session_state.show_back_toast = True
            #set game mode to 'single_setup' to return to single player setup page
            st.session_state.game_mode = "single_setup"
            #force a rerun to apply the change immediately
            st.rerun()


if __name__ == "__main__":
    st.set_page_config(
        page_title="Molecular Pictionary",
        layout="centered",
        initial_sidebar_state="collapsed",
    )
    render_game_page()