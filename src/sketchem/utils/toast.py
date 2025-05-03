#This is an implementation of toast notifications for streamlit, similar to how they work in something like Next js

import streamlit as st
import time

def create_toast(message, type="info", duration=5):
    """
    Create a toast notification that persists across reruns
    
    Args:
        message: The message to display
        type: The type of toast (info, success, error, warning)
        duration: How long the toast should remain visible in seconds
    """
    # Store the toast in session state with expiration time
    current_time = time.time()
    st.session_state.toast = {
        "message": message,
        "type": type,
        "expires_at": current_time + duration
    }

def show_toast():
    """
    Display any active toast notifications and remove expired ones

    Should be called at the beginning of each page for the toast to display above everything else
    """
    if "toast" in st.session_state:
        toast = st.session_state.toast
        st.session_state.current_time = time.time()
        
        # Check if toast is still valid
        if st.session_state.current_time < toast["expires_at"]:
            # Display the toast based on its type
            if toast["type"] == "info":
                st.info(toast["message"])
            elif toast["type"] == "success":
                st.success(toast["message"])
            elif toast["type"] == "error":
                st.error(toast["message"])
            elif toast["type"] == "warning":
                st.warning(toast["message"])
                
            # Add custom styling to make it look more like a toast
            st.markdown("""
            <style>
            div[data-testid="stInfoMessage"],
            div[data-testid="stSuccessMessage"],
            div[data-testid="stWarningMessage"],
            div[data-testid="stErrorMessage"] {
                position: fixed;
                bottom: 20px;
                right: 20px;
                z-index: 9999;
                min-width: 300px;
                max-width: 400px;
                animation: fadeIn 0.5s;
            }
            @keyframes fadeIn {
                from { opacity: 0; }
                to { opacity: 1; }
            }
            </style>
            """, unsafe_allow_html=True)
        else:
            # Toast has expired, remove it
            del st.session_state.toast
            st.rerun() #Rerun the page to remove the toast??