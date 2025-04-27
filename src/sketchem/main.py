import streamlit as st
from .utils.game_state import initialize_game_state
from .pages.home_page import render_home_page



def main():
    st.set_page_config(page_title="Sketchem", layout="centered")
    st.title("🧪 Sketchem")

    # Initialize the game state parameters
    initialize_game_state()

    # Route to appropriate page based on game mode chosen
    if st.session_state.game_mode is None:
        render_home_page()
    elif st.session_state.game_mode == "singleplayer_setup":
        # Render setup for single here
        pass
    elif st.session_state.game_mode == "multiplayer_setup":
        # Render setup for multiplayer here
        pass
    elif st.session_state.game_mode == "single":
        # Render single player game here
        pass
    elif st.session_state.game_mode == "multiplayer":
        # Render multiplayer game here
        pass

if __name__ == "__main__":
    main()