import streamlit as st
from streamlit_drawable_canvas import st_canvas
from PIL import Image
import io
from streamlit_extras.vertical_slider import vertical_slider


def save_canvas_as_image(canvas_data):  # convert canvas data to png image
    if canvas_data is not None:
        img_data = canvas_data.astype("uint8")
        img = Image.fromarray(img_data[..., :3])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None


# This function was replaced by direct pen/eraser buttons


def render_game_page():
    # initialize session states
    if "pen_size" not in st.session_state:
        st.session_state.pen_size = 3
    if "eraser_size" not in st.session_state:
        st.session_state.eraser_size = 10
    if "canvas_key" not in st.session_state:
        st.session_state.canvas_key = 0
    if "drawing_mode" not in st.session_state:
        st.session_state.drawing_mode = "freedraw"  # default to pen mode
    if "last_pen_color" not in st.session_state:
        st.session_state.last_pen_color = "White"  # default
    if "device_mode" not in st.session_state:
        st.session_state.device_mode = "desktop"  # default

    # color options
    color_options = {
        "White":  "#ffffff",
        "Red":    "#ff0000",
        "Blue":   "#0000ff",
        "Green":  "#00ff00",
        "Yellow": "#ffff00",
    }

    # configure canvas based on mode
    if st.session_state.drawing_mode == "erase":
        current_stroke_color = "#000000"
        current_stroke_width = st.session_state.eraser_size if "eraser_size" in st.session_state else 10
    else:
        current_stroke_color = color_options[st.session_state.last_pen_color]
        current_stroke_width = st.session_state.pen_size

    #title
    st.markdown(
        """
    <h1 style='text-align: center; margin-bottom: 40px; color: #ffffff;'>Single Player Mode</h1>
    <div style='margin-top: 40px;'></div>
    """,
        unsafe_allow_html=True,
    )

    # Device selection
    device_cols = st.columns([1, 1])
    with device_cols[0]:
        if st.button("💻", key="desktop_btn", use_container_width=True):
            st.session_state.device_mode = "desktop"
            st.rerun()

    with device_cols[1]:
        if st.button("📱", key="mobile_btn", use_container_width=True):
            st.session_state.device_mode = "mobile"
            st.rerun()

    # Highlight active device button
    st.markdown(f"""
    <style>
    div[data-testid="stButton"] button[key="desktop_btn"] {{
        {("border: 2px solid #4CAF50 !important;") if st.session_state.device_mode == "desktop" else ""}
    }}
    div[data-testid="stButton"] button[key="mobile_btn"] {{
        {("border: 2px solid #4CAF50 !important;") if st.session_state.device_mode == "mobile" else ""}
    }}
    </style>
    """, unsafe_allow_html=True)

    is_mobile = st.session_state.device_mode == "mobile"

    # Add CSS to align controls with canvas
    st.markdown("""
    <style>
    /* Align controls with canvas */
    [data-testid="column"]:first-child {
        align-self: flex-start;
        padding-top: 0;
    }
    </style>
    """, unsafe_allow_html=True)

    #layout: controls vs. canvas
    col1, col2 = st.columns([1, 3])

    #controls
    with col1:
        # Create two columns for pen and eraser buttons
        tool_cols = st.columns(2)

        # Pen button
        with tool_cols[0]:
            if st.button("✏️", key="pen_button"):
                if st.session_state.drawing_mode == "erase":
                    st.session_state.drawing_mode = "freedraw"
                    st.rerun()

        # Eraser button
        with tool_cols[1]:
            if st.button("🧽", key="eraser_button"):
                if st.session_state.drawing_mode == "freedraw":
                    st.session_state.drawing_mode = "erase"
                    st.rerun()

        # Active button highlight
        st.markdown("""
        <style>
        div[data-testid="stButton"] button[key="pen_button"] {
            %s
        }
        div[data-testid="stButton"] button[key="eraser_button"] {
            %s
        }
        </style>
        """ % (
            "border: 2px solid #4CAF50 !important;" if st.session_state.drawing_mode == "freedraw" else "",
            "border: 2px solid #ff6b6b !important;" if st.session_state.drawing_mode == "erase" else ""
        ), unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        if st.session_state.drawing_mode != "erase":

            #color buttons
            color_cols = st.columns(len(color_options))

            # Create colored buttons
            for i, (color_name, color_code) in enumerate(color_options.items()):
                with color_cols[i]:
                    # Create a button for each color
                    if st.button("",
                               key=f"color_btn_{color_name}",
                               help=color_name,
                               use_container_width=True):
                        # Update the pen color immediately when clicked
                        st.session_state.last_pen_color = color_name
                        st.rerun()

            # Add CSS to style the buttons with their respective colors
            for color_name, color_code in color_options.items():
                is_selected = color_name == st.session_state.last_pen_color
                border = "3px solid #4CAF50" if is_selected else "2px solid #888"

                st.markdown(f"""
                <style>
                /* Style for {color_name} button */
                div[data-testid="stButton"] button[key="color_btn_{color_name}"] {{
                    background-color: {color_code} !important;
                    border: {border} !important;
                    border-radius: 5px !important;
                    height: 30px !important;
                    padding: 0 !important;
                }}
                </style>
                """, unsafe_allow_html=True)

            # Update current stroke color
            current_stroke_color = color_options[st.session_state.last_pen_color]

        st.markdown("<br>", unsafe_allow_html=True)

        #size slider: different for desktop and mobile
        if st.session_state.drawing_mode == "erase":
            # Eraser size slider
            eraser_size = st.slider(
                "",
                min_value=5,
                max_value=30,
                value=st.session_state.eraser_size,
                key="eraser_size_slider"
            )

            #update immediately if changed
            if eraser_size != st.session_state.eraser_size:
                st.session_state.eraser_size = eraser_size
                st.rerun()

            current_stroke_width = eraser_size
        else:
            #pen size slider: vertical for desktop, horizontal for mobile
            if is_mobile:
                pen_size = st.slider(
                    "",
                    min_value=1,
                    max_value=20,
                    value=st.session_state.pen_size,
                    key="pen_size_slider"
                )

                #update immediately if changed
                if pen_size != st.session_state.pen_size:
                    st.session_state.pen_size = pen_size
                    st.rerun()
            else:
                #vertical slider for desktop
                try:
                    # Get the canvas height to match the slider height
                    canvas_height = 100  # Same as the canvas height

                    pen_size = vertical_slider(
                        label="",
                        min_value=1,
                        max_value=20,
                        default_value=st.session_state.pen_size,
                        key="pen_size_slider",
                        height=canvas_height  # Match the canvas height
                    )

                    #udate immediately if changed
                    if pen_size != st.session_state.pen_size:
                        st.session_state.pen_size = pen_size
                        st.rerun()
                except Exception:
                    #fallback to regular slider
                    pen_size = st.slider(
                        "",
                        min_value=1,
                        max_value=20,
                        value=st.session_state.pen_size,
                        key="fallback_pen_size_slider"
                    )

                    #update immediately if changed
                    if pen_size != st.session_state.pen_size:
                        st.session_state.pen_size = pen_size
                        st.rerun()

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
            drawing_mode=st.session_state.drawing_mode,
            key=f"canvas_{st.session_state.canvas_key}",
            display_toolbar=True,
        )

    #buttons row - add more space below the canvas
    st.markdown("<div style='margin-top: 80px;'></div>", unsafe_allow_html=True)

    # 2 equal columns for buttons
    submit_col, back_col = st.columns([1, 1])

    #submit button (left)
    with submit_col:
        if st.button("Submit Drawing", type="primary", key="submit_btn", use_container_width=True):
            if canvas_result.image_data is not None:
                img_bytes = save_canvas_as_image(canvas_result.image_data)
                if img_bytes:
                    st.toast("Drawing submitted successfully!", icon="✅")
            else:
                st.toast("Please draw something before submitting!", icon="⚠️")

    #back button (right)
    with back_col:
        # Use a simple button instead of the back_button function
        if st.button("Back", key="simple_back_btn", use_container_width=True):
            st.session_state.show_back_toast = True
            st.session_state.game_mode = "single_setup"  # Go back to single player setup
            st.rerun()


if __name__ == "__main__":
    st.set_page_config(
        page_title="Molecular Pictionary",
        layout="centered",
        initial_sidebar_state="collapsed",
    )
    render_game_page()