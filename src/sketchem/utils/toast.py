"""
Toast notification utilities for Sketchem.

This file contains functions to display toast notifications
that persist across page reloads.
"""

import streamlit as st

def display_queued_toast():
    if st.session_state.get("toast_queue"):
        message = st.session_state.toast_queue.get("message", "No message")
        icon = st.session_state.toast_queue.get("icon") #st.toast handles None icon okay
        st.toast(message, icon=icon)
        #clear the queue immediately after displaying
        st.session_state.toast_queue = None
