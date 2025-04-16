from fastapi import FastAPI, WebSocket
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi import Request
import uvicorn
import logging
import random
from .game_manager import GameManager

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("sketchem_app")

app = FastAPI()
logger.info("Initializing GameManager")
game_manager = GameManager()

# Mount static files
logger.info("Mounting static files")
app.mount("/static", StaticFiles(directory="src/static"), name="static")
templates = Jinja2Templates(directory="src/templates")
logger.info("App configuration complete")

@app.get("/", response_class=HTMLResponse)
async def get(request: Request):
    logger.info("Serving main page")
    return templates.TemplateResponse("layout.html", {"request": request})

@app.get("/reset", response_class=RedirectResponse)
async def reset(request: Request):
    logger.info("Reset requested - clearing game manager state")
    # Use the reset method instead of recreating
    game_manager.reset()
    # Redirect to home with a cache-busting parameter
    return RedirectResponse(url=f"/?nocache={random.randint(1000, 9999)}")

@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    client_host = websocket.client.host
    logger.info(f"New WebSocket connection from {client_host}")
    try:
        await game_manager.handle_connection(websocket)
    except Exception as e:
        logger.error(f"Error handling WebSocket connection: {e}", exc_info=True)

if __name__ == "__main__":
    logger.info("Starting Sketchem server")
    uvicorn.run("src.main:app", host="0.0.0.0", port=8000, reload=True) 