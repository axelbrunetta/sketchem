# This is used to install the DECIMER model if it is not already installed

import subprocess
import sys
from pathlib import Path
import os
import tempfile
from .environment import is_running_locally
import streamlit as st

def get_decimer_path():
    """Get the appropriate DECIMER path based on environment"""
    if is_running_locally():
        return Path.home() / '.decimer'  # Local path
    # Use a temporary directory for Streamlit Cloud
    return Path(tempfile.gettempdir()) / 'decimer'

@st.cache_resource
def install_decimer_model():
    """Install the DECIMER Canonical model if not already installed"""
    try:
        # Use environment-specific path
        decimer_path = get_decimer_path() / 'DECIMER_model_weights' / 'Canonical'
        
        if decimer_path.exists():
            return True
            
        # Ensure directory exists
        os.makedirs(get_decimer_path(), exist_ok=True)
            
        print("Installing DECIMER Canonical model...")
        result = subprocess.run(
            [sys.executable, "-m", "decimer", "--model", "Canonical", "--image", "dummy"],
            capture_output=True,
            text=True,
            env=dict(os.environ, DECIMER_DATA_PATH=str(get_decimer_path()))
        )
        
        if result.returncode != 0:
            print(f"Error installing DECIMER model: {result.stderr}")
            return False
            
        return True
    except Exception as e:
        print(f"Error during DECIMER model installation: {str(e)}")
        return False
