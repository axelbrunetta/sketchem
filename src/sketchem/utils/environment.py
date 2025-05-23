"""
Environment utilities for Sketchem.

This file contains helper functions for detecting the runtime environment
and retrieving configuration values like API keys.
"""

import platform
import os
import streamlit as st



def is_running_locally():
    """
    Check if the app is running locally or in Streamlit cloud.
    
    Returns:
        Boolean indicating if the app is running locally
    """
    return platform.processor() != ''

def get_gemini_api_key():
    """
    Get the Google Gemini API key from environment variables or Streamlit secrets.
    
    Tries to get the key from Streamlit secrets first, then falls back to
    loading from .env file if running locally.
    
    Returns:
        The Gemini API key as a string
    """
    try:
        if "GEMINI_API_KEY" in st.secrets:
            return st.secrets["GEMINI_API_KEY"]
    except Exception:
        from dotenv import load_dotenv
        load_dotenv()  # Load environment variables from .env file
        return os.environ.get("GEMINI_API_KEY")
