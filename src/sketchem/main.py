import streamlit as st
from sketchem.utils.game_state import initialize_game_state
from sketchem.pages.home_page import render_home_page
from sketchem.pages.single_player_setup import render_singleplayer_setup
from sketchem.pages.multiplayer_setup import render_multiplayer_setup
from sketchem.pages.waiting_room import render_waiting_room
from sketchem.utils.toast import display_queued_toast
from sketchem.pages.game_page_single import render_game_page
from sketchem.pages.game_page_multi import render_game_page_multi
from sketchem.pages.guide_page import render_guide_page
from streamlit_js_eval import streamlit_js_eval

def main():
    st.set_page_config(page_title="Sketchem", layout="wide", initial_sidebar_state="collapsed")

    # Get actual screen width -> used for mobile-specific styling
    st.session_state.viewport_width = streamlit_js_eval(js_expressions="window.innerWidth", key="test_viewport_width")

    st.session_state.is_mobile = st.session_state.viewport_width < 768 if st.session_state.viewport_width else False
    
    st.markdown(
        """
    <style>
        [data-testid="collapsedControl"] {
            display: none
        }
        section[data-testid="stSidebar"] {
            display: none;
        }
        div[data-testid="stSidebarCollapsedControl"] {
            display: none;
        }
    </style>
    """,
        unsafe_allow_html=True,
    )
    display_queued_toast() #Show any active toast notifications

    col1, col2, col3 = st.columns([1, 2, 1]) # Columns are now needed because "wide" mode has to be enabled for streamlit so that the proper width of the screen can be determined for phone / computer detection -> otherwise everything is stretched by wide mode
    

    #with col2:
        #st.title("sketchem")

    # Initialize the game state parameters
    initialize_game_state()

    # Route to appropriate page based on game mode chosen
    if st.session_state.game_mode is None:
        render_home_page()
    elif st.session_state.game_mode == "singleplayer_setup":
        render_singleplayer_setup()
    elif st.session_state.game_mode == "multiplayer_setup": # reroute to multiplayer setup page
        render_multiplayer_setup()
    elif st.session_state.game_mode in ["created_multi", "joined_multi"]: # reroute to waiting room for both host and joining players
        render_waiting_room()
    elif st.session_state.game_mode == "single":
        render_game_page()
    elif st.session_state.game_mode == "multiplayer":
        render_game_page_multi()
    elif st.session_state.game_mode == "guide":
        render_guide_page()


if __name__ == "__main__":
    main()
