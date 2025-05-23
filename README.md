![Project Logo](assets/banner.png)

<h1 align="center">
Sketchem
</h1>

<br>


Sketchem is an interactive molecular pictionary game where users compete to draw molecules correctly. The app uses AI to recognize hand-drawn molecular structures and provides a fun way to learn about chemistry through gameplay.

## Table of Contents

- [Live Game](#Live-game)
- [Features](#Features)
- [Installation](#Installation)
- [Requirements](#Requirements)
- [How to Run Locally](#How-to-Run-Locally)
- [How to deploy on Streamlit Cloud](#How-to-deploy-on-Streamlit-Cloud)
- [Run Tests and Coverage](#Run-Tests-and-Coverage)
- [Troubleshooting](#Troubleshooting)
- [Contributions](#Contributions)
- [License](#License)

## Live Game

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://sketchem.streamlit.app)

## Features

- **Single-player and multiplayer modes**: Challenge yourself or compete with friends in real time
- **Paper-like drawing canvas**: Sketch molecules with an intuitive drawing interface
- **AI-powered recognition**: Google Generative AI (Gemini) interprets hand-drawn structures to identify the drawn molecule
- **AI-powered molecule category generation**: Gemini AI creates categories for molecules based on user input
- **Leaderboard system**: Track scores based on accuracy and speed, and compete for the top spot!
- **Hints**: Get feedback on your drawings to improve your structure-drawing skills

## Installation

```bash
# Create a new environment
conda create -n sketchem python=3.10
conda activate sketchem

# Clone the repository
git clone https://github.com/axelbrunetta/sketchem.git
cd sketchem

# Install the package
pip install .
```


## Requirements


### Gemini AI

To use the Gemini AI for molecule recognition and category creation, you need to set up your Google Generative AI API key:

1. Go to [Google Generative AI](https://aistudio.google.com/app/apikey) and create an API key
2. Create a `.env` file in the root of your project with:
   ```
   GEMINI_API_KEY=your-gemini-api-key
   ```

### Jupyter (Optional)

If you need Jupyter Lab to open the notebook:

```bash
(sketchem) $ pip install jupyterlab
```

## How to Run Locally

1. Make sure your `.env` file is in the root of your project with (from requirements above):
   ```
   GEMINI_API_KEY=your-gemini-api-key
   ```

2. Run the app:
   ```bash
   streamlit run src/sketchem/main.py
   ```

**Note**: Multiplayer game mode is not available when running the app locally; it only works on the deployed version.

## How to deploy on Streamlit Cloud


1. Fork or push this repository to your GitHub account
2. Connect your repository to Streamlit Cloud
3. Add the following secret in the Streamlit Cloud dashboard:
```
GEMINI_API_KEY = "your-gemini-api-key"
```
4. In you app's settings, change the python version to < 3.12 as 3.12+ leads to issues with dependencies



### Run Tests

```bash
(sketchem) $ pip install tox
(sketchem) $ tox
```



## Troubleshooting

### StreamlitSetPageConfigMustBeFirstCommandError

If you encounter this error when deploying on Streamlit Cloud:
```
streamlit.errors.StreamlitSetPageConfigMustBeFirstCommandError: This app has encountered an error.
```

This is a normal occurrence during the first load of the app on Streamlit Cloud and is linked to us disabling the default Streamlit sidebar. Simply reload the page once and the error will disappear. 

## Contributions


Ivana and Ariadna worked on the singleplayer game page, the guide page, and styling of various pages.

Axel worked on the multiplayer integration, the multiplayer game page, and other features like category creation and molecule recognition.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
