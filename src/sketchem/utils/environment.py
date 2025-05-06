"""Environment utilities for Sketchem."""

import platform
import os
import streamlit as st

# Behavior for streamlit is different in cloud and locally run -> especially for what we will do with the database, hence the need for this

def is_running_locally():
    """
    Determine if the application is running locally or in a deployed environment.

    Returns:
        bool: True if running locally, False if running in a deployed environment
    """
    # This works since Streamlit Cloud runs on Linux without a processor name
    return platform.processor() != ''

def get_gemini_api_key():
    """Get Gemini API key from environment"""
    try:
        if "GEMINI_API_KEY" in st.secrets:
            return st.secrets["GEMINI_API_KEY"]
    except Exception:
        return os.environ.get("GEMINI_API_KEY")
