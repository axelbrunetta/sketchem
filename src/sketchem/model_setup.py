# This file is used to install the DECIMER model if it is not already installed

import subprocess
import sys
from pathlib import Path

def install_decimer_model():
    """Install the DECIMER Canonical model if not already installed."""
    try:
        # Check if model is already installed
        decimer_path = Path.home() / '.decimer' / 'DECIMER_model_weights' / 'Canonical'
        if decimer_path.exists():
            return True
            
        print("Installing DECIMER Canonical model...")
        result = subprocess.run(
            [sys.executable, "-m", "decimer", "--model", "Canonical", "--image", "dummy"],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"Error installing DECIMER model: {result.stderr}")
            return False
            
        return True
    except Exception as e:
        print(f"Error during DECIMER model installation: {str(e)}")
        return False