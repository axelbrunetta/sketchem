"""Run the Sketchem application."""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.abspath("src"))

# Import and run the main function
from sketchem.main import main

if __name__ == "__main__":
    main()
