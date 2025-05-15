"""Environment utilities for Sketchem."""

import platform
import os
import streamlit as st
from dotenv import load_dotenv
from pathlib import Path

# Behavior for streamlit is different in cloud and locally run -> especially for what we will do with the database, hence the need for this


def is_running_locally():
    """
    Determine if the application is running locally or in a deployed environment.

    Returns:
        bool: True if running locally, False if running in a deployed environment
    """
    # This works since Streamlit Cloud runs on Linux without a processor name
    return platform.processor() != ""


def get_gemini_api_key():
    """Get Gemini API key from environment (.env locally or Streamlit secrets in deployed)."""
    # Local environment: load from .env file
    if is_running_locally():
        # Determine project root and .env path
        root_dir = Path(__file__).resolve().parents[3]
        env_path = root_dir / ".env"
        if not env_path.exists():
            raise ValueError(
                f".env file not found at {env_path}. Please create it and add your GEMINI_API_KEY."
            )
        load_dotenv(dotenv_path=env_path, override=True)
        api_key = os.getenv("GEMINI_API_KEY")
        if not api_key:
            raise ValueError(
                "GEMINI_API_KEY not found in .env file. Please add it to the .env file in the root directory."
            )

        return api_key.strip()
    # Deployed environment: load from Streamlit secrets or OS environment
    try:
        if "GEMINI_API_KEY" in st.secrets:
            return st.secrets["GEMINI_API_KEY"]
    except Exception:
        pass
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        raise ValueError(
            "GEMINI_API_KEY not found in Streamlit secrets or environment variables."
        )
    return api_key.strip()