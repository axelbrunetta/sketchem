"""Mock database for Sketchem."""

import uuid
import random
import string
import time
from streamlit.logger import get_logger
import logging

logger = get_logger(__name__)
logger.setLevel(logging.DEBUG)

# In-memory database for games
GAMES = {}

def generate_game_code():
    """Generate a random 4-character game code."""
    return ''.join(random.choices(string.ascii_uppercase, k=4))

def create_game(player_name):
    """Create a new game and return the game code."""
    game_code = generate_game_code()
    
    # Ensure the game code is unique
    while game_code in GAMES:
        game_code = generate_game_code()
    
    # Generate a unique player ID
    player_id = str(uuid.uuid4())
    
    # Create the game
    GAMES[game_code] = {
        "created_at": time.time(),
        "started": False,
        "players": {
            player_id: {
                "name": player_name,
                "score": 0,
                "is_host": True
            }
        }
    }
    
    logger.info(f"Created game {game_code} with host {player_name}")
    
    return {
        "game_code": game_code,
        "player_id": player_id
    }

def join_game(game_code, player_name):
    """Join an existing game."""
    # Check if the game exists
    if game_code not in GAMES:
        return {"success": False, "error": "Game not found"}
    
    # Check if the game has already started
    if GAMES[game_code]["started"]:
        return {"success": False, "error": "Game has already started"}
    
    # Generate a unique player ID
    player_id = str(uuid.uuid4())
    
    # Add the player to the game
    GAMES[game_code]["players"][player_id] = {
        "name": player_name,
        "score": 0,
        "is_host": False
    }
    
    logger.info(f"Player {player_name} joined game {game_code}")
    
    return {
        "success": True,
        "player_id": player_id
    }

def get_game(game_code):
    """Get a game by its code."""
    return GAMES.get(game_code)

def start_game(game_code):
    """Start a game."""
    # Check if the game exists
    if game_code not in GAMES:
        return {"success": False, "error": "Game not found"}
    
    # Mark the game as started
    GAMES[game_code]["started"] = True
    
    logger.info(f"Game {game_code} started")
    
    return {"success": True}
