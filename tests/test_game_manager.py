import pytest
from src.game_manager import GameManager
from fastapi import WebSocket

@pytest.mark.asyncio
async def test_join_room():
    game_manager = GameManager()
    websocket = WebSocket()
    room_id = "test_room"
    
    await game_manager.join_room(websocket, room_id)
    
    assert room_id in game_manager.rooms
    assert id(websocket) in game_manager.rooms[room_id]
    assert id(websocket) in game_manager.active_connections

@pytest.mark.asyncio
async def test_disconnect():
    game_manager = GameManager()
    websocket = WebSocket()
    room_id = "test_room"
    
    await game_manager.join_room(websocket, room_id)
    await game_manager.disconnect(websocket)
    
    assert room_id not in game_manager.rooms
    assert id(websocket) not in game_manager.active_connections

@pytest.mark.asyncio
async def test_broadcast_draw():
    game_manager = GameManager()
    websocket1 = WebSocket()
    websocket2 = WebSocket()
    room_id = "test_room"
    
    await game_manager.join_room(websocket1, room_id)
    await game_manager.join_room(websocket2, room_id)
    
    message = {
        "type": "draw",
        "data": {
            "x1": 0,
            "y1": 0,
            "x2": 100,
            "y2": 100
        }
    }
    
    # Mock the send_json method
    websocket2.send_json = pytest.AsyncMock()
    
    await game_manager.broadcast_draw(websocket1, message)
    
    # Verify that websocket2 received the message
    websocket2.send_json.assert_called_once_with(message) 