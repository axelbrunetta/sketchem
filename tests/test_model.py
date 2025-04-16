import pytest
from src.ml_model import MLModel
import numpy as np

def test_model_initialization():
    model = MLModel()
    assert model.model is None
    assert not model.is_loaded

def test_model_loading():
    model = MLModel()
    model.load_model("test_model_path")
    assert model.is_loaded

def test_preprocess_drawing():
    model = MLModel()
    drawing_data = {
        "strokes": [
            {"x": [0, 100], "y": [0, 100]},
            {"x": [100, 200], "y": [100, 200]}
        ]
    }
    processed_data = model.preprocess_drawing(drawing_data)
    assert isinstance(processed_data, np.ndarray)

def test_predict_without_loading():
    model = MLModel()
    with pytest.raises(RuntimeError):
        model.predict({"test": "data"})

def test_model_info():
    model = MLModel()
    info = model.get_model_info()
    assert isinstance(info, dict)
    assert "is_loaded" in info
    assert "model_type" in info
    assert "version" in info 