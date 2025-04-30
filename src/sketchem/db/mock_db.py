import uuid
import time
import streamlit as st
from typing import Dict, List, Optional
import random
import string
from streamlit.logger import get_logger
import logging

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

# In-memory storage / Mock database -> Could replace this with an actual database like Firestore but not needed here since Streamlit cloud only runs one instance of the code hence everyone can use the "same storage"

_games = {} # This will be a dictionary of all the games running

def generate_game_code(length: int = 6) -> str: #Self-explanatory
    """Generate a random game code"""
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length)) 

def create_game(player_name: str) -> Dict:
    """Create a new game with a unique code"""

    code = generate_game_code() #Generate code using function above

    while code in _games: # Check that game does not already exist / isn't overwritten
        code = generate_game_code()

    player_id = str(uuid.uuid4()) # Generate unique code (UUID) for each player in case players choose to have the same name -> would create conflict since this we're using common storage
    
    game_data = {
        "code": code,
        "status": "waiting",
        "created_at": int(time.time()),
        "category": st.session_state.selected_molecule_category,
        "categoryIsDefault": st.session_state.categoryIsDefault,
        "additionalCategories": st.session_state.additionalCategories, #Adds additional categories to database when creating the game
        "game_duration": st.session_state.game_duration,
        "players": {
            player_id: {
                "name": player_name,
                "id": player_id,
                "score": 0,
                "last_active": int(time.time()),
            }
        },
    }
    
    _games[code] = game_data # Add game element to the fake database using the code as the key and the game data as the attached value 
    return {"game_code": code, "player_id": player_id}



def join_game(code: str, player_name: str) -> Dict: #Processes code entered by player and allowsthem to join a hosted game 
    """Join an existing game"""
    if code not in _games:
        logger.info("Game not found")
        return {"success": False, "error": "Game not found"}
        
    game = _games[code]
    logger.info(f"Game: {game}")

    player_id = str(uuid.uuid4()) #No need to check that player name alr exists as we're using unique identifiers
    game["players"][player_id] = { 
        "name": player_name,
        "id": player_id,
        "score": 0,
        "last_active": int(time.time())
    }
    logger.info("Successfully joined game")
    return {"success": True, "player_id": player_id}

def get_game(code: str) -> Optional[Dict]: #Returns game for given code -> used in waiting room to display game info
    """Get game data"""
    if code not in _games:
        return None
    return _games[code]



def start_game(code: str) -> Dict: #Self-explanatory
    """Start the game"""
    if code not in _games:
        return {"success": False, "error": "Game not found"}
        
    _games[code]["status"] = "active"
    return {"success": True}

#Will need end game function here?

#Will need something that deletes games and players after x amount of time