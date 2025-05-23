
import platform
import os
import streamlit as st


# Behavior for streamlit is different in cloud and locally run -> especially for what we will do with the databse, hence the need for this

def is_running_locally():
    """Check if the app is running locally or in Streamlit cloud -> this works since Streamlit Cloud runs on Linux without a processor name"""
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
        load_dotenv()  #load environment variables from .env file
        return os.environ.get("GEMINI_API_KEY")
