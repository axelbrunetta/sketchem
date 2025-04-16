# Sketchem: Molecule Drawing Game

A multiplayer web-based game where players compete to draw molecules correctly.

## Features

- Create and join game rooms
- Real-time drawing with WebSocket communication
- ML-based verification of molecule drawings
- Multiplayer scoreboard and round tracking
- Responsive design that works on desktop and mobile

## Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/sketchem.git
cd sketchem
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
python -m src.main
```

4. Open a web browser and navigate to `http://localhost:8000`

## How to Play

1. Create a new game room or join an existing one with a room code
2. When all players have joined, the host can start the game
3. Each round, players will see a molecule to draw
4. Draw the molecule on the canvas and submit your drawing
5. The ML model will check if your drawing is correct
6. Score points for correct drawings
7. The player with the highest score after all rounds wins!

## Technical Details

- Backend: FastAPI + WebSockets
- Frontend: HTML5 Canvas + JavaScript
- ML: Uses a pre-trained model to identify molecule drawings

## Adding Your Own Molecules

To add your own molecules to the game:

1. Add molecule images to the `src/static/molecules/` directory
2. Update the list of molecules in `src/game_manager.py`

## Testing

Run the tests using pytest:

```bash
pytest tests/
``` 