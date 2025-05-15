import streamlit as st
from sketchem.utils.environment import is_running_locally
from streamlit_extras.bottom_container import bottom
from streamlit_extras.stoggle import stoggle


def render_home_page():
    st.markdown("## Choose Game Mode")

    col1, col2 = st.columns(2)

    with col1:
        if st.button("Single Player", use_container_width=True):
            st.session_state.game_mode = "singleplayer_setup"
            st.rerun()

    with col2:
        if (
            not is_running_locally()
        ):  # Prevents multiplayer button from showing when running locally
            if st.button("Multiplayer", use_container_width=True):
                st.session_state.game_mode = "multiplayer_setup"
                st.rerun()
        else:
            st.info(
                "Multiplayer is only available in the deployed version (Using Streamlit Cloud)"
            )

    if st.button("How to play", use_container_width=True):
        st.session_state.game_mode = "guide"
        st.rerun()

    with bottom():
        stoggle(
            "Credits",
            """Created by Axel Brunetta, Ivana Josipovic and Ariadna Davila""",
        )
