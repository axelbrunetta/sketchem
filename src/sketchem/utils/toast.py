#This is an implementation of toast notifications for streamlit, similar to how they work in something like Next js / React

import streamlit as st
import time
import threading
from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx
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
    
    st.session_state.toast = {
        "message": message,
        "type": type,
    }
    
    # Create a thread to handle the delayed toast removal
    class ToastRemovalThread(threading.Thread):
        def __init__(self, delay):
            super().__init__()
            self.delay = delay
            
        def run(self):
            time.sleep(self.delay)
            if "toast" in st.session_state:
                del st.session_state.toast
                st.rerun()
    
    # Get current context and create thread with context
    ctx = get_script_run_ctx()
    if ctx is not None:
        thread = ToastRemovalThread(duration)
        add_script_run_ctx(thread, ctx)
        thread.start()
    
def show_toast():
    """
    Display any active toast notifications and remove expired ones

    Should be called at the beginning of each page for the toast to display above everything else
    """
    # Check if toast exists in session state
    if "toast" in st.session_state:
        toast = st.session_state.toast

        # Display the toast based on its type
        if toast["type"] == "info":
            st.info(toast["message"])
        elif toast["type"] == "success":
            st.success(toast["message"])
        elif toast["type"] == "error":
            st.error(toast["message"])
        elif toast["type"] == "warning":
            st.warning(toast["message"])
                
            
        
