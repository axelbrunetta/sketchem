"""A molecular pictionary."""

from __future__ import annotations

__version__ = "0.0.1"

from .model_setup import install_decimer_model

# Install DECIMER model on first import
install_decimer_model()

# Import your package's functions
from .smiles_validator import validate_drawing