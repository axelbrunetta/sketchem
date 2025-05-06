# This is used to install the DECIMER model if it is not already installed

import subprocess
import sys
from pathlib import Path
import os
from .environment import is_running_locally

def get_decimer_path():
    """Get the appropriate DECIMER path based on environment"""
    if is_running_locally():
        return Path.home() / '.decimer' # Local path
    return Path('/mount/decimer')  # Streamlit Cloud path

def install_decimer_model():
    """Install the DECIMER Canonical model if not already installed"""
    try:
        # Use environment-specific path
        decimer_path = get_decimer_path() / 'DECIMER_model_weights' / 'Canonical'

        if decimer_path.exists():
            print("DECIMER model already installed.")
            return True

        if not is_running_locally():
            # In cloud, ensure mount directory exists
            os.makedirs('/mount/decimer', exist_ok=True)

        print("Installing DECIMER Canonical model...")

        # First try to import DECIMER to check if it's installed
        try:
            import DECIMER
            print("DECIMER module found, attempting to download model...")

            # Try to use the DECIMER API directly instead of command line
            try:
                # This should trigger the model download
                DECIMER.predict_SMILES("dummy", "Canonical")
                print("DECIMER model installed successfully.")
                return True
            except Exception as model_error:
                print(f"Error using DECIMER API: {str(model_error)}")
                # Fall back to command line approach
        except ImportError:
            print("DECIMER module not found in Python path.")

        # Fall back to command line approach
        try:
            result = subprocess.run(
                [sys.executable, "-m", "decimer", "--model", "Canonical", "--image", "dummy"],
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                print(f"Error installing DECIMER model via command line: {result.stderr}")
                return False

            return True
        except Exception as cmd_error:
            print(f"Error running DECIMER command: {str(cmd_error)}")
            return False

    except Exception as e:
        print(f"Error during DECIMER model installation: {str(e)}")
        return False
