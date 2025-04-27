# This file contains the tests for smile_validator.py

import pytest
from pathlib import Path
from sketchem.utils.smiles_validator import validate_drawing

# Get current working directory and construct test data path
current_directory = Path(__file__).parent.parent
TEST_DATA_DIR = current_directory / 'data' / 'testdata'

def test_morgan_comparison():
    """Test Morgan fingerprint comparison with different thresholds"""
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
    """Test handling of invalid inputs"""
    image_path = TEST_DATA_DIR / 'ethanol.png'
    
    # Test invalid method
    with pytest.raises(ValueError):
        validate_drawing(str(image_path), 'CCO', method='invalid')
    
    # Test invalid SMILES
    assert validate_drawing(str(image_path), 'invalid_smiles') == False
    
    # Test invalid image path
    assert validate_drawing('nonexistent.png', 'CCO') == False

def test_threshold_bounds():
    """Test Morgan comparison with edge case thresholds"""
    image_path = TEST_DATA_DIR / 'ethanol.png'
    
    # Test threshold = 1.0 (exact match required)
    result_exact = validate_drawing(str(image_path), 'CCO', 
                                  method='morgan', threshold=1.0)
    assert result_exact == True
    
    # Test very low threshold
    result_low = validate_drawing(str(image_path), 'CCO', 
                                method='morgan', threshold=0.1)
    assert result_low == True

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

if __name__ == '__main__':
    pytest.main([__file__])
