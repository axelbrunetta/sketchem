import streamlit as st
from sketchem.utils.game_state import initialize_game_state
from sketchem.pages.home_page import render_home_page
from sketchem.pages.multiplayer_setup import render_multiplayer_setup
from sketchem.pages.waiting_room import render_waiting_room
from sketchem.utils.toast import display_queued_toast
from sketchem.pages.game_page_single import render_game_page
from sketchem.pages.game_page_multi import render_game_page_multi

def main():
    st.set_page_config(page_title="Sketchem", layout="centered", initial_sidebar_state="collapsed")

    st.markdown(
        """
    <style>
        [data-testid="collapsedControl"] {
            display: none
        }
        section[data-testid="stSidebar"] {
            display: none;
        }
    </style>
    """,
        unsafe_allow_html=True,
    )
    display_queued_toast() #Show any active toast notifications

    st.title("🧪 Sketchem")

    # Initialize the game state parameters
    initialize_game_state()

    # Route to appropriate page based on game mode chosen
    if st.session_state.game_mode is None:
        render_home_page()
    elif st.session_state.game_mode == "singleplayer_setup":
        # Render setup for single here
        pass
    elif st.session_state.game_mode == "multiplayer_setup": # reroute to multiplayer setup page
        render_multiplayer_setup()
    elif st.session_state.game_mode in ["created_multi", "joined_multi"]: # reroute to waiting room for both host and joining players
        render_waiting_room()
    elif st.session_state.game_mode == "single":
        render_game_page()
    elif st.session_state.game_mode == "multiplayer":
        render_game_page_multi()

if __name__ == "__main__":
    main()
