
import platform
import os
import streamlit as st

# Behavior for streamlit is different in cloud and locally run -> especially for what we will do with the databse, hence the need for this

def is_running_locally():
    """Check if the app is running locally or in Streamlit cloud -> this works since Streamlit Cloud runs on Linux without a processor name"""
    return platform.processor() != ''

def get_gemini_api_key():
    """Get Gemini API key from environment"""
    try:
        if "GEMINI_API_KEY" in st.secrets:
            return st.secrets["GEMINI_API_KEY"]
    except Exception:
        return os.environ.get("GEMINI_API_KEY")