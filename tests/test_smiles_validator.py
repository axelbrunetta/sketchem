# This file contains the tests for smile_validator.py

import pytest
from pathlib import Path
from sketchem.utils.smiles_validator import validate_drawing, validate_drawing_with_ai
from sketchem.utils.environment import get_gemini_api_key

# Get current working directory and construct test data path
current_directory = Path(__file__).parent.parent
TEST_DATA_DIR = current_directory / 'data' / 'testdata'

from dotenv import load_dotenv
load_dotenv()

"""
def test_morgan_comparison():
    
    # Test cases with expected results
    test_cases = [
        {
            'image': 'ethanol.png',
            'target': 'CCO',  # ethanol
            'similar': 'CCCO',  # propanol (similar structure)
            'different': 'c1ccccc1',  # benzene (different structure)
        },
        # Add more test cases as needed
    ]
    
    for case in test_cases:
        image_path = TEST_DATA_DIR / case['image']
        
        # Test with different thresholds
        # High threshold (0.85) - strict matching
        result_strict = validate_drawing(str(image_path), case['target'], 
                                      method='morgan', threshold=0.85)
        assert result_strict == True
        
        # Medium threshold (0.75) - should match similar structures
        result_medium = validate_drawing(str(image_path), case['similar'], 
                                       method='morgan', threshold=0.75)
        assert result_medium == True
        
        # Different molecule should fail even with low threshold
        result_different = validate_drawing(str(image_path), case['different'], 
                                          method='morgan', threshold=0.65)
        assert result_different == False


def test_invalid_inputs():
    
    image_path = TEST_DATA_DIR / 'ethanol.png'
    
    # Test invalid method
    with pytest.raises(ValueError):
        validate_drawing(str(image_path), 'CCO', method='invalid')
    
    # Test invalid SMILES
    assert validate_drawing(str(image_path), 'invalid_smiles') == False
    
    # Test invalid image path
    assert validate_drawing('nonexistent.png', 'CCO') == False

def test_threshold_bounds():
    
    image_path = TEST_DATA_DIR / 'ethanol.png'
    
    # Test threshold = 1.0 (exact match required)
    result_exact = validate_drawing(str(image_path), 'CCO', 
                                  method='morgan', threshold=1.0)
    assert result_exact == True
    
    # Test very low threshold
    result_low = validate_drawing(str(image_path), 'CCO', 
                                method='morgan', threshold=0.1)
    assert result_low == True 
"""

def test_mcs_comparison():
    """Test Maximum Common Substructure comparison with different thresholds"""
    # Test cases with expected results
    test_cases = [
        {
            'image': 'ethanol.png',
            'target': 'CCO',  # ethanol
            'similar': 'CCCO',  # propanol (similar structure)
            'different': 'c1ccccc1',  # benzene (different structure)
        },
        # Add more test cases as needed
    ]
    
    for case in test_cases:
        image_path = TEST_DATA_DIR / case['image']
        
        # Test with different thresholds
        # High threshold (0.8) - strict matching
        result_strict = validate_drawing(str(image_path), case['target'], 
                                      method='mcs', threshold=0.8)
        assert result_strict == True, "MCS strict matching failed"
        
        # Medium threshold (0.7) - should match similar structures
        result_medium = validate_drawing(str(image_path), case['similar'], 
                                       method='mcs', threshold=0.7)
        assert result_medium == True, "MCS medium threshold matching failed"
        
        # Different molecule should fail even with low threshold
        result_different = validate_drawing(str(image_path), case['different'], 
                                          method='mcs', threshold=0.6)
        assert result_different == False, "MCS different structure check failed"

def test_mcs_edge_cases():
    """Test MCS comparison with edge cases"""
    image_path = TEST_DATA_DIR / 'ethanol.png'
    
    # Test with maximum threshold
    result_max = validate_drawing(str(image_path), 'CCO', 
                                method='mcs', threshold=1.0)
    assert result_max == True, "MCS maximum threshold failed"
    
    # Test with very low threshold
    result_low = validate_drawing(str(image_path), 'CCO', 
                                method='mcs', threshold=0.1)
    assert result_low == True, "MCS minimum threshold failed"
    
    # Test with invalid SMILES
    result_invalid = validate_drawing(str(image_path), 'invalid_smiles', 
                                    method='mcs')
    assert result_invalid == False, "MCS invalid SMILES handling failed"

def test_validate_drawing_with_ai():
    """Test AI-based validation using Gemini"""

    image_path = TEST_DATA_DIR / 'ethanol.png'
    api_key = get_gemini_api_key()
    
    if not api_key:
        pytest.skip("Gemini API key not found")
    
    # Test cases with expected results
    test_cases = [
        {
            'target': 'CCO',  # ethanol -> should match
            'threshold': 0.85,
            'expected': True
        },
        {
            'target': 'c1ccccc1',  # benzene -> should not match
            'threshold': 0.85,
            'expected': False
        }
    ]

    for case in test_cases:
        result = validate_drawing_with_ai(
            api_key=api_key,
            image_path=str(image_path),
            target_smiles=case['target'],
            threshold=case['threshold']
        )
        assert result == case['expected'], f"Expected {case['expected']} for {case['target']}"
        assert isinstance(result, bool), "Result should be boolean" # We also verify the function runs without errors and returns expected type

def test_validate_drawing_with_ai_invalid_inputs():
    """Test AI validation with invalid inputs"""

    image_path = TEST_DATA_DIR / 'ethanol.png'
    api_key = get_gemini_api_key()
    
    # Test with invalid API key
    result = validate_drawing_with_ai(
        api_key="",
        image_path=str(image_path),
        target_smiles='CCO'
    )
    assert result == "❗ Gemini API key not set."
    
    if api_key:  # Only run these tests if we have an API key
        # Test with invalid SMILES
        result = validate_drawing_with_ai(
            api_key=api_key,
            image_path=str(image_path),
            target_smiles='invalid_smiles'
        )
        assert result is False

        # Test with invalid image path - only test if we have API key
        result = validate_drawing_with_ai(
            api_key=api_key,
            image_path='nonexistent.png',
            target_smiles='CCO'
        )
        assert isinstance(result, str) and "error" in result.lower()

if __name__ == '__main__':
    pytest.main([__file__])
