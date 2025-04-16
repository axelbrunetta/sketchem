class MoleculeDrawingGame {
    constructor() {
        console.log("Initializing MoleculeDrawingGame v1.2 - Fixed Game Over Panel");
        
        // Game elements
        this.setupPanel = document.getElementById('game-setup');
        this.roomPanel = document.getElementById('game-room');
        this.resultOverlay = document.getElementById('result-overlay');
        this.gameOverPanel = document.getElementById('game-over');
        
        // Force hiding game over panel at construction time
        this.forceHideGameOver();
        
        // Canvas elements
        this.canvas = document.getElementById('drawing-canvas');
        this.ctx = this.canvas.getContext('2d');
        
        // Game state
        this.isDrawing = false;
        this.lastX = 0;
        this.lastY = 0;
        this.ws = null;
        this.roomId = null;
        this.playerId = null;
        this.gameState = null;
        this.isHost = false;
        this.reconnectAttempts = 0;
        this.maxReconnectAttempts = 3;

        this.setupCanvas();
        this.setupEventListeners();
        
        // Double-check that game over is hidden after everything is initialized
        setTimeout(() => this.forceHideGameOver(), 100);
    }
    
    forceHideGameOver() {
        // Force game over panel to be hidden with multiple techniques
        if (this.gameOverPanel) {
            console.log("Forcing game over panel to be hidden");
            this.gameOverPanel.classList.add('hidden');
            this.gameOverPanel.style.display = 'none';
            document.querySelector('#game-over').classList.add('hidden');
            document.querySelector('#game-over').style.display = 'none';
        } else {
            console.error("Could not find game-over element");
        }
    }

    setupCanvas() {
        // Make the canvas responsive based on its container
        this.resizeCanvas();
        window.addEventListener('resize', () => this.resizeCanvas());
        
        // Set drawing style
        this.ctx.strokeStyle = '#000';
        this.ctx.lineWidth = 2;
        this.ctx.lineCap = 'round';
        this.ctx.lineJoin = 'round';
        
        // Clear any existing content
        this.clearCanvas();
        
        // Ensure canvas is visible
        const canvasContainer = document.querySelector('.canvas-container');
        canvasContainer.classList.remove('hidden');
        this.canvas.style.display = 'block';
    }
    
    resizeCanvas() {
        const container = this.canvas.parentElement;
        const rect = container.getBoundingClientRect();
        this.canvas.width = rect.width;
        this.canvas.height = rect.height - document.querySelector('.tools').clientHeight;
        console.log('Canvas resized to:', this.canvas.width, 'x', this.canvas.height);
    }

    setupEventListeners() {
        // Room buttons
        document.getElementById('create-room').addEventListener('click', () => this.createRoom());
        document.getElementById('join-room').addEventListener('click', () => this.joinRoom());
        document.getElementById('start-game').addEventListener('click', () => this.startGame());
        
        // Canvas drawing events
        this.canvas.addEventListener('mousedown', this.startDrawing.bind(this));
        this.canvas.addEventListener('mousemove', this.draw.bind(this));
        this.canvas.addEventListener('mouseup', this.stopDrawing.bind(this));
        this.canvas.addEventListener('mouseout', this.stopDrawing.bind(this));

        // Touch events for mobile
        this.canvas.addEventListener('touchstart', this.handleTouchStart.bind(this));
        this.canvas.addEventListener('touchmove', this.handleTouchMove.bind(this));
        this.canvas.addEventListener('touchend', this.stopDrawing.bind(this));

        // Game buttons
        document.getElementById('clear-canvas').addEventListener('click', () => this.clearCanvas());
        document.getElementById('submit-drawing').addEventListener('click', () => this.submitDrawing());
        document.getElementById('continue-game').addEventListener('click', () => this.hideResultOverlay());
        document.getElementById('new-game').addEventListener('click', () => this.resetGame());
    }

    createRoom() {
        this.connectWebSocket(() => {
            this.ws.send(JSON.stringify({
                type: 'create'
            }));
        });
    }

    joinRoom() {
        const roomId = document.getElementById('room-id').value;
        if (!roomId) {
            alert('Please enter a room code');
            return;
        }

        this.connectWebSocket(() => {
            this.ws.send(JSON.stringify({
                type: 'join',
                room_id: roomId
            }));
        });
    }

    startGame() {
        if (this.ws && this.ws.readyState === WebSocket.OPEN) {
            this.ws.send(JSON.stringify({
                type: 'start_game'
            }));
        }
    }

    connectWebSocket(onOpenCallback) {
        console.log('Attempting WebSocket connection...');
        const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
        const wsUrl = `${protocol}//${window.location.host}/ws`;
        
        this.ws = new WebSocket(wsUrl);

        this.ws.onopen = () => {
            console.log('WebSocket connection established');
            this.reconnectAttempts = 0; // Reset reconnect attempts on successful connection
            if (onOpenCallback) onOpenCallback();
        };

        this.ws.onmessage = (event) => {
            try {
                const message = JSON.parse(event.data);
                console.log('Received WebSocket message:', message);
                this.handleServerMessage(message);
            } catch (e) {
                console.error('Error parsing WebSocket message:', e);
            }
        };

        this.ws.onerror = (error) => {
            console.error('WebSocket error:', error);
        };

        this.ws.onclose = (event) => {
            console.log('WebSocket connection closed:', event.code, event.reason);
            
            // Ensure game over panel is hidden when connection closes
            if (this.gameOverPanel) {
                this.gameOverPanel.classList.add('hidden');
            }
            
            // Only attempt to reconnect if we have a room ID (meaning we were in a game)
            if (this.roomId && this.reconnectAttempts < this.maxReconnectAttempts) {
                this.reconnectAttempts++;
                console.log(`Attempting to reconnect (${this.reconnectAttempts}/${this.maxReconnectAttempts})...`);
                setTimeout(() => this.reconnect(), 1000 * this.reconnectAttempts);
            }
        };
    }

    reconnect() {
        console.log('Attempting to reconnect to room:', this.roomId);
        if (this.isHost) {
            this.connectWebSocket(() => {
                this.ws.send(JSON.stringify({
                    type: 'create',
                    room_id: this.roomId
                }));
            });
        } else {
            this.connectWebSocket(() => {
                this.ws.send(JSON.stringify({
                    type: 'join',
                    room_id: this.roomId
                }));
            });
        }
    }

    handleServerMessage(message) {
        console.log('Received message:', message);
        
        switch (message.type) {
            case 'room_created':
                console.log('Handling room_created message:', message);
                this.handleRoomCreated(message);
                break;
            case 'room_joined':
                console.log('Handling room_joined message:', message);
                this.handleRoomJoined(message);
                break;
            case 'room_state':
                console.log('Handling room_state message:', message);
                this.handleRoomState(message.state);
                break;
            case 'error':
                console.error('Received error:', message.message);
                alert(message.message);
                break;
            case 'prediction_result':
                console.log('Handling prediction_result:', message);
                this.handlePredictionResult(message);
                break;
            case 'draw':
                console.log('Handling draw message:', message);
                this.handleRemoteDrawing(message.data);
                break;
        }
    }

    handleRoomCreated(message) {
        console.log('Room created with ID:', message.room_id);
        this.roomId = message.room_id;
        this.playerId = message.player_id;
        this.isHost = true;
        
        // Force hide game over panel
        this.forceHideGameOver();
        
        // Initialize game state for new room
        this.gameState = {
            players: {[this.playerId]: {name: "Player 1", score: 0}},
            current_molecule: null,
            round: 0,
            max_rounds: 5,
            status: "waiting"
        };
        console.log('Initialized game state:', this.gameState);
        
        // Switch to room panel
        this.setupPanel.classList.add('hidden');
        this.roomPanel.classList.remove('hidden');
        
        // Update room code display
        document.getElementById('room-code').textContent = this.roomId;
        
        // Show start button for host
        document.getElementById('start-game').classList.remove('hidden');
        
        // Update display with initial state
        this.updateGameDisplay();
    }

    handleRoomJoined(message) {
        this.roomId = message.room_id;
        this.playerId = message.player_id;
        this.gameState = message.room_state;
        
        // Switch to room panel
        this.setupPanel.classList.add('hidden');
        this.roomPanel.classList.remove('hidden');
        this.gameOverPanel.classList.add('hidden'); // Ensure game over panel is hidden
        this.resultOverlay.classList.add('hidden'); // Ensure result overlay is hidden
        
        // Update room code display
        document.getElementById('room-code').textContent = this.roomId;
        
        // Update player list
        this.updateScoreboard();
        
        // Update game display
        this.updateGameDisplay();
    }

    handleRoomState(state) {
        console.log('Updating room state:', state);
        this.gameState = state;
        this.updateGameDisplay();
    }
    
    handlePredictionResult(result) {
        const resultMessage = document.getElementById('result-message');
        
        if (result.correct) {
            resultMessage.textContent = 'Correct! Moving to next molecule.';
            resultMessage.style.color = '#2ecc71';
        } else {
            resultMessage.textContent = `Incorrect. Try again. The molecule was ${result.expected}.`;
            resultMessage.style.color = '#e74c3c';
        }
        
        this.resultOverlay.classList.remove('hidden');
    }

    updateGameDisplay() {
        if (!this.gameState) {
            console.log('No game state available');
            return;
        }
        
        console.log('Updating game display with state:', this.gameState);
        
        // Update status
        const statusElement = document.getElementById('room-status');
        const statusText = this.gameState.status === 'lobby' || this.gameState.status === 'waiting'
            ? 'Waiting for host to start game...' 
            : this.gameState.status === 'playing' 
                ? 'Game in progress' 
                : 'Game finished';
        statusElement.textContent = statusText;
        
        // Update scoreboard and round info
        this.updateScoreboard();
        document.getElementById('current-round').textContent = this.gameState.round;
        document.getElementById('max-rounds').textContent = this.gameState.max_rounds;
        
        // Show/hide start button for host
        const startButton = document.getElementById('start-game');
        if (this.isHost && (this.gameState.status === 'lobby' || this.gameState.status === 'waiting')) {
            startButton.classList.remove('hidden');
        } else {
            startButton.classList.add('hidden');
        }
        
        // Update current molecule
        this.updateMoleculeDisplay();
        
        // Handle canvas visibility
        const canvasContainer = document.querySelector('.canvas-container');
        const lobbyMessage = document.getElementById('lobby-message') || this.createLobbyMessage();
        
        if (this.gameState.status === 'lobby' || this.gameState.status === 'waiting') {
            canvasContainer.classList.add('hidden');
            lobbyMessage.classList.remove('hidden');
            
            const message = this.isHost 
                ? 'You are the host. Click "Start Game" when everyone has joined.'
                : 'Waiting for the host to start the game...';
            lobbyMessage.querySelector('p').textContent = message;
        } else {
            canvasContainer.classList.remove('hidden');
            lobbyMessage.classList.add('hidden');
            this.setupCanvas(); // Reinitialize canvas when game starts
        }
    }
    
    updateMoleculeDisplay() {
        const molecule = this.gameState.current_molecule;
        const nameElement = document.getElementById('molecule-name');
        const formulaElement = document.getElementById('molecule-formula');
        const imageElement = document.getElementById('molecule-image');
        
        if (molecule && this.gameState.status === 'playing') {
            nameElement.textContent = molecule.name;
            formulaElement.textContent = molecule.formula;
            imageElement.src = molecule.image_url;
            imageElement.alt = `${molecule.name} structure`;
        } else {
            nameElement.textContent = 'Waiting for game to start...';
            formulaElement.textContent = '';
            imageElement.src = '';
            imageElement.alt = '';
        }
    }
    
    updateScoreboard() {
        if (!this.gameState || !this.gameState.players) return;
        
        const playersListElement = document.getElementById('players-list');
        playersListElement.innerHTML = '';
        
        Object.entries(this.gameState.players).forEach(([id, player]) => {
            const playerElement = document.createElement('div');
            playerElement.className = 'player-score';
            if (id === this.playerId) {
                playerElement.classList.add('current-player');
            }
            
            playerElement.innerHTML = `
                <span class="player-name">${player.name}</span>
                <span class="score">${player.score}</span>
            `;
            
            playersListElement.appendChild(playerElement);
        });
    }
    
    showGameOver() {
        console.log('Preparing to show game over panel');
        
        const finalScoresElement = document.getElementById('final-scores');
        finalScoresElement.innerHTML = '<h3>Final Scores</h3>';
        
        // Sort players by score
        const sortedPlayers = Object.entries(this.gameState.players)
            .sort(([, a], [, b]) => b.score - a.score);
        
        sortedPlayers.forEach(([id, player], index) => {
            const playerElement = document.createElement('div');
            playerElement.className = 'player-score';
            if (id === this.playerId) {
                playerElement.classList.add('current-player');
            }
            
            playerElement.innerHTML = `
                <span class="player-name">${index + 1}. ${player.name}</span>
                <span class="score">${player.score}</span>
            `;
            
            finalScoresElement.appendChild(playerElement);
        });
        
        // Only now, after everything is prepared, show the game over panel
        this.gameOverPanel.classList.remove('hidden');
        this.gameOverPanel.style.display = 'flex';
        
        console.log('Game over panel should now be visible');
    }
    
    hideResultOverlay() {
        this.resultOverlay.classList.add('hidden');
        this.clearCanvas();
    }
    
    resetGame() {
        // Hide game over panel
        this.gameOverPanel.classList.add('hidden');
        
        // Go back to setup screen
        this.roomPanel.classList.add('hidden');
        this.setupPanel.classList.remove('hidden');
        
        // Reset game state
        this.roomId = null;
        this.playerId = null;
        this.gameState = null;
        this.isHost = false;
        
        // Close WebSocket connection
        if (this.ws) {
            this.ws.close();
            this.ws = null;
        }
    }

    getCoordinates(e) {
        const rect = this.canvas.getBoundingClientRect();
        const scaleX = this.canvas.width / rect.width;
        const scaleY = this.canvas.height / rect.height;
        
        // Get the actual mouse position relative to the canvas
        const x = (e.clientX - rect.left) * scaleX;
        const y = (e.clientY - rect.top) * scaleY;
        
        console.log('Mouse position:', e.clientX, e.clientY);
        console.log('Canvas rect:', rect);
        console.log('Calculated coordinates:', x, y);
        
        return [x, y];
    }

    startDrawing(e) {
        if (this.gameState?.status !== 'playing') {
            console.log('Cannot draw: game not in playing state');
            return;
        }
        
        this.isDrawing = true;
        const [x, y] = this.getCoordinates(e);
        this.lastX = x;
        this.lastY = y;
        console.log('Started drawing at:', x, y);
    }

    draw(e) {
        if (!this.isDrawing) return;

        const [x, y] = this.getCoordinates(e);
        console.log('Drawing to:', x, y);
        
        this.ctx.beginPath();
        this.ctx.moveTo(this.lastX, this.lastY);
        this.ctx.lineTo(x, y);
        this.ctx.stroke();

        // Send drawing data to server
        if (this.ws && this.ws.readyState === WebSocket.OPEN) {
            this.ws.send(JSON.stringify({
                type: 'draw',
                data: {
                    x1: this.lastX,
                    y1: this.lastY,
                    x2: x,
                    y2: y
                }
            }));
        }
        
        [this.lastX, this.lastY] = [x, y];
    }
    
    handleRemoteDrawing(data) {
        this.ctx.beginPath();
        this.ctx.moveTo(data.x1, data.y1);
        this.ctx.lineTo(data.x2, data.y2);
        this.ctx.stroke();
    }

    stopDrawing() {
        this.isDrawing = false;
    }

    handleTouchStart(e) {
        e.preventDefault();
        this.startDrawing(e.touches[0]);
    }

    handleTouchMove(e) {
        e.preventDefault();
        this.draw(e.touches[0]);
    }

    clearCanvas() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    }

    submitDrawing() {
        if (!this.ws || this.ws.readyState !== WebSocket.OPEN || this.gameState?.status !== 'playing') {
            return;
        }
        
        // Get the image data from canvas
        const imageData = this.canvas.toDataURL('image/png');
        
        // Send to server for verification
        this.ws.send(JSON.stringify({
            type: 'predict',
            data: {
                image: imageData,
                molecule: this.gameState.current_molecule.name
            }
        }));
    }

    // Create a lobby message element if it doesn't exist
    createLobbyMessage() {
        const container = document.createElement('div');
        container.id = 'lobby-message';
        container.className = 'lobby-message';
        
        const message = document.createElement('p');
        message.textContent = 'Waiting for the game to start...';
        container.appendChild(message);
        
        // Add after canvas container
        const canvasContainer = document.querySelector('.canvas-container');
        canvasContainer.parentNode.insertBefore(container, canvasContainer.nextSibling);
        
        return container;
    }
}

// Initialize the game when the page loads
window.addEventListener('load', () => {
    new MoleculeDrawingGame();
}); 