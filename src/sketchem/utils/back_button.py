import streamlit as st
from sketchem.utils.game_state import reset_game_state
from streamlit_extras import stylable_container

def back_button(destination=None, label="Back", use_container_width=True, key=None):
    """
    Renders a back button that navigates to the specified destination or previous page.
    
    Args:
        destination: The game_mode to return to. If None, will try to determine the previous page.
        label: in case you want to change the label to "Home" or something else
        use_container_width: Whether the button should use the full container width
        key: Optional key for the button
    
    Returns:
        True if the button was clicked, False otherwise
    """
    button_key = key or f"back_button_{destination}" # buttons and elements in streamlit need unique identifier -> prevents errors if multiple back buttons because of a problem with the refresh
    
    # Create columns with 0.1 and 0.9 width ratio
    
    
    
    with stylable_container(
    key="Upload_Data",
    css_styles="""
    button{
        float: right;
    }
    """ ):
        if st.button(label, use_container_width=use_container_width, key=button_key):
            if destination is not None:
                # Navigate to specific destination
                st.session_state.game_mode = destination
            else:
                # Determine previous page based on current game_mode
                current_mode = st.session_state.game_mode
                
                # Define reroutes
                navigation_map = { #for now all back buttons go to home page for simplicity
                    "multiplayer": None,  # multiplayer game goes back to home
                    "multiplayer_setup": None,  # setup goes back to home
                    "created_multi": None,  # waiting room (host) goes back home
                    "joined_multi": None,  # waiting room (player) goes back home
                    # need one for multiplayer goes home from game page (exits game)
                    # need one for single player goes back to home from setup + game page
                }
                
                # Set the previous page as the destination
                st.session_state.game_mode = navigation_map.get(current_mode, None)
            
            # Clear any page-specific state if needed
            if current_mode in ["created_multi", "joined_multi", "multiplayer"]:
                reset_game_state()
            
            st.rerun()
            return True
    
    # Return false if button wasn't clicked
    return False
