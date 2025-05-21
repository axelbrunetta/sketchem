import streamlit as st
from sketchem.utils.environment import is_running_locally
from streamlit_extras.bottom_container import bottom
from streamlit_extras.stoggle import stoggle
from streamlit_js_eval import streamlit_js_eval
import os

                   
def render_home_page():
    #get logo path
    current_dir = os.path.dirname(__file__)
    logo_path = os.path.join(current_dir, "..", "..", "..", "assets", "logo.jpg")

    #create columns
    padding1, goodcolumn, padding2 = st.columns([1, 2, 1])
    
    #display logo
    with goodcolumn:
        image_width = 300
        st.image(logo_path, width=image_width if not st.session_state.is_mobile else None, use_container_width=True if st.session_state.is_mobile else False)

        #display title and buttons  
        st.markdown("## Choose Game Mode")
        col4, col5 = st.columns(2)
        with col4:
            if st.button("Single Player", use_container_width=True):
                st.session_state.game_mode = "singleplayer_setup"
                st.rerun()

        with col5:
            if not is_running_locally(): # Prevents multiplayer button from showing when running locally
                if st.button("Multiplayer", use_container_width=True):
                    st.session_state.game_mode = "multiplayer_setup"
                    st.rerun()
            else:
                st.info("Multiplayer is only available in the deployed version (Using Streamlit Cloud)")

        if st.button("How to play", use_container_width=True):
                    st.session_state.game_mode = "guide"
                    st.rerun()

    #display credits
    with bottom():
        stoggle(
        "Credits",
        """Created by Axel Brunetta, Ivana Josipovic and Ariadna Davila""",
        )