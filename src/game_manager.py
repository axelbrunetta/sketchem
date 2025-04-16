from typing import Dict, Set, List, Optional
from fastapi import WebSocket
import json
import random
from .ml_model import MLModel

class GameManager:
    def __init__(self):
        self.active_connections: Dict[str, WebSocket] = {}
        self.rooms: Dict[str, Set[str]] = {}
        self.room_states: Dict[str, Dict] = {}
        self.ml_model = MLModel()
        self.ml_model.load_model("model_path")  # You'll need to specify your model path
        self.molecules = [
            {"name": "Water", "formula": "H2O", "image_url": "/static/molecules/water.png"},
            {"name": "Carbon Dioxide", "formula": "CO2", "image_url": "/static/molecules/co2.png"},
            {"name": "Methane", "formula": "CH4", "image_url": "/static/molecules/methane.png"},
            {"name": "Glucose", "formula": "C6H12O6", "image_url": "/static/molecules/glucose.png"},
            {"name": "Benzene", "formula": "C6H6", "image_url": "/static/molecules/benzene.png"}
        ]
        print("GameManager initialized with empty state")

    async def handle_connection(self, websocket: WebSocket):
        print("New WebSocket connection established")
        await websocket.accept()
        try:
            while True:
                data = await websocket.receive_text()
                message = json.loads(data)
                await self.handle_message(websocket, message)
        except Exception as e:
            print(f"Error handling connection: {e}")
        finally:
            await self.disconnect(websocket)

    async def handle_message(self, websocket: WebSocket, message: dict):
        print(f"Received message: {message}")
        message_type = message.get("type")
        if message_type == "join":
            room_id = message.get("room_id")
            if room_id:
                print(f"Joining room {room_id}")
                await self.join_room(websocket, room_id)
        elif message_type == "create":
            print("Creating new room...")
            await self.create_room(websocket)
        elif message_type == "draw":
            print("Broadcasting draw event")
            await self.broadcast_draw(websocket, message)
        elif message_type == "predict":
            print("Verifying drawing")
            await self.verify_drawing(websocket, message.get("data"))
        elif message_type == "start_game":
            print("Starting game")
            await self.start_game(websocket)

    async def create_room(self, websocket: WebSocket):
        print("Starting room creation...")
        # Generate a unique room ID
        room_id = f"room_{random.randint(1000, 9999)}"
        while room_id in self.rooms:
            room_id = f"room_{random.randint(1000, 9999)}"
        
        print(f"Generated room ID: {room_id}")
        
        # Initialize room
        self.rooms[room_id] = set()
        connection_id = str(id(websocket))
        self.rooms[room_id].add(connection_id)
        self.active_connections[connection_id] = websocket
        
        # Initialize room state
        self.room_states[room_id] = {
            "players": {connection_id: {"name": "Player 1", "score": 0}},
            "current_molecule": None,
            "round": 0,
            "max_rounds": 5,
            "status": "waiting"
        }
        
        print(f"Room state initialized: {self.room_states[room_id]}")
        
        # Send room info back to creator
        await websocket.send_json({
            "type": "room_created",
            "room_id": room_id,
            "player_id": connection_id
        })
        print("Room created message sent to client")

    async def join_room(self, websocket: WebSocket, room_id: str):
        if room_id not in self.rooms:
            await websocket.send_json({
                "type": "error",
                "message": "Room does not exist"
            })
            return
            
        if self.room_states[room_id]["status"] != "waiting":
            await websocket.send_json({
                "type": "error",
                "message": "Game already in progress"
            })
            return
            
        connection_id = str(id(websocket))
        self.rooms[room_id].add(connection_id)
        self.active_connections[connection_id] = websocket
        
        # Update player count
        player_number = len(self.room_states[room_id]["players"]) + 1
        self.room_states[room_id]["players"][connection_id] = {
            "name": f"Player {player_number}",
            "score": 0
        }
        
        # Send room info to the new player
        await websocket.send_json({
            "type": "room_joined",
            "room_id": room_id,
            "player_id": connection_id,
            "room_state": self.room_states[room_id]
        })
        
        # Notify all players about the new player
        await self.broadcast_room_state(room_id)

    async def start_game(self, websocket: WebSocket):
        connection_id = str(id(websocket))
        room_id = self.get_room_for_connection(connection_id)
        
        if not room_id:
            print(f"Cannot start game: connection {connection_id} not in any room")
            return
            
        if self.room_states[room_id]["status"] != "waiting":
            print(f"Cannot start game: room {room_id} is not in 'waiting' status")
            return
            
        print(f"Starting game in room {room_id}")
        # Start the game with the first molecule
        self.room_states[room_id]["status"] = "playing"
        self.room_states[room_id]["round"] = 1
        self.room_states[room_id]["current_molecule"] = random.choice(self.molecules)
        
        print(f"Room state after starting game: {self.room_states[room_id]}")
        
        # Broadcast game start and first molecule
        await self.broadcast_room_state(room_id)

    async def verify_drawing(self, websocket: WebSocket, drawing_data):
        connection_id = str(id(websocket))
        room_id = self.get_room_for_connection(connection_id)
        
        if not room_id:
            print(f"Cannot verify drawing: connection {connection_id} not in any room")
            return
            
        if self.room_states[room_id]["status"] != "playing":
            print(f"Cannot verify drawing: room {room_id} is not in 'playing' status")
            return
            
        current_molecule = self.room_states[room_id]["current_molecule"]
        if not current_molecule:
            print(f"Cannot verify drawing: no current molecule in room {room_id}")
            return
            
        print(f"Verifying drawing for molecule {current_molecule['name']} in room {room_id}")
        
        # Use ML model to verify if the drawing matches the molecule
        try:
            # For debugging: pretend the prediction matches to see if game flow works correctly
            # Comment this out and uncomment the line below when ML model is implemented
            is_correct = True  # For testing
            # is_correct = self.ml_model.predict(drawing_data) == current_molecule["name"]
            
            print(f"Prediction result: {'correct' if is_correct else 'incorrect'}")
            
            # Update player score if correct
            if is_correct:
                self.room_states[room_id]["players"][connection_id]["score"] += 10
                print(f"Player {connection_id} score increased to {self.room_states[room_id]['players'][connection_id]['score']}")
                
                # Move to next round
                self.room_states[room_id]["round"] += 1
                print(f"Moving to round {self.room_states[room_id]['round']}")
                
                if self.room_states[room_id]["round"] > self.room_states[room_id]["max_rounds"]:
                    # Game over
                    print(f"Game over in room {room_id} - reached max rounds")
                    self.room_states[room_id]["status"] = "finished"
                    self.room_states[room_id]["current_molecule"] = None
                else:
                    # Select next molecule
                    self.room_states[room_id]["current_molecule"] = random.choice(self.molecules)
                    print(f"Selected new molecule: {self.room_states[room_id]['current_molecule']['name']}")
        except Exception as e:
            print(f"Error in verify_drawing: {e}")
            is_correct = False
        
        # Send result to the player
        await websocket.send_json({
            "type": "prediction_result",
            "correct": is_correct,
            "expected": current_molecule["name"]
        })
        
        # Broadcast updated state to all players
        await self.broadcast_room_state(room_id)

    async def broadcast_draw(self, websocket: WebSocket, message: dict):
        connection_id = str(id(websocket))
        room_id = self.get_room_for_connection(connection_id)
        
        if room_id:
            for conn_id in self.rooms[room_id]:
                if conn_id != connection_id:
                    await self.active_connections[conn_id].send_json(message)

    async def broadcast_room_state(self, room_id: str):
        for connection_id in self.rooms[room_id]:
            await self.active_connections[connection_id].send_json({
                "type": "room_state",
                "state": self.room_states[room_id]
            })

    def get_room_for_connection(self, connection_id: str) -> Optional[str]:
        for room_id, connections in self.rooms.items():
            if connection_id in connections:
                return room_id
        return None

    async def disconnect(self, websocket: WebSocket):
        # Remove from active connections
        connection_id = str(id(websocket))
        if connection_id in self.active_connections:
            del self.active_connections[connection_id]
        
        # Remove from rooms and update room states
        room_to_delete = None
        for room_id, connections in self.rooms.items():
            if connection_id in connections:
                connections.remove(connection_id)
                
                # Remove player from room state
                if room_id in self.room_states and connection_id in self.room_states[room_id]["players"]:
                    del self.room_states[room_id]["players"][connection_id]
                
                # If room is empty, mark for deletion
                if not connections:
                    room_to_delete = room_id
                else:
                    # Otherwise broadcast updated state
                    await self.broadcast_room_state(room_id)
        
        # Delete empty room
        if room_to_delete:
            del self.rooms[room_to_delete]
            if room_to_delete in self.room_states:
                del self.room_states[room_to_delete]

    def reset(self):
        """Reset all game state"""
        print("Resetting all game state")
        self.active_connections = {}
        self.rooms = {}
        self.room_states = {}
        print("Game state has been reset") 