import streamlit as st
from sketchem.utils.game_state import initialize_game_state
from sketchem.pages.home_page import render_home_page
from sketchem.pages.multiplayer_setup import render_multiplayer_setup
from sketchem.pages.waiting_room import render_waiting_room
from sketchem.pages.single_player_setup import render_singleplayer_setup
from sketchem.pages.game_page_single import render_game_page as render_single_game


def main():
    st.set_page_config(
        page_title="Sketchem",
        layout="centered",
        menu_items={},
        initial_sidebar_state="collapsed"
    )

    # Initialize the game state parameters
    initialize_game_state()

    # Route to appropriate page based on game mode chosen
    if st.session_state.game_mode is None:
        render_home_page()
    elif st.session_state.game_mode == "single_setup":
        # Render setup for single player
        render_singleplayer_setup()
    elif st.session_state.game_mode == "multiplayer_setup":
        render_multiplayer_setup()
    elif st.session_state.game_mode in ["created_multi", "joined_multi"]:
        render_waiting_room()
    elif st.session_state.game_mode == "single":
        # Render single player game
        render_single_game()
    elif st.session_state.game_mode == "multiplayer":
        # Render multiplayer game here
        pass

if __name__ == "__main__":
    main()
