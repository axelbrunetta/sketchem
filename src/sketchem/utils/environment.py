"""Environment utilities for Sketchem."""

import os
import streamlit as st

def is_running_locally():
    """
    Determine if the application is running locally or in a deployed environment.
    
    Returns:
        bool: True if running locally, False if running in a deployed environment
    """
    # For our purposes, we'll always return False to enable the multiplayer button
    return False
