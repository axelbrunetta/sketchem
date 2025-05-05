# Setting Up Gemini API for Sketchem

This guide will help you set up the Gemini API key for the Sketchem application.

## Getting a Gemini API Key

1. Go to the [Google AI Studio](https://makersuite.google.com/app/apikey) website.
2. Sign in with your Google account.
3. Click on "Get API key" or "Create API key".
4. Copy the generated API key.

## Adding the API Key to Sketchem

1. In your Sketchem project directory, create a `.streamlit` folder if it doesn't exist:
   ```bash
   mkdir -p .streamlit
   ```

2. Create a `secrets.toml` file inside the `.streamlit` folder:
   ```bash
   touch .streamlit/secrets.toml
   ```

3. Open the `secrets.toml` file and add your Gemini API key:
   ```toml
   # .streamlit/secrets.toml
   GEMINI_API_KEY = "your-api-key-here"
   ```
   Replace `"your-api-key-here"` with the actual API key you obtained from Google AI Studio.

4. Save the file and restart the Streamlit application.

## Testing the API Connection

1. Run the Sketchem application:
   ```bash
   streamlit run src/sketchem/main.py
   ```

2. Navigate to the Single Player or Multiplayer setup page.
3. Try creating a custom molecule category using the "Create a molecule category" button.
4. If the API key is set up correctly, you should be able to generate custom molecule categories.

## Troubleshooting

- If you see an error message saying "Gemini API key not found in secrets", check that your `secrets.toml` file is in the correct location and has the correct format.
- Make sure the API key is valid and has not expired.
- Check that you have the required Python packages installed:
  ```bash
  pip install google-generativeai
  ```

## Security Notes

- Never commit your `secrets.toml` file to version control.
- The `.streamlit` folder is already included in the `.gitignore` file to prevent accidental commits.
- If you're deploying the application, use the appropriate method for your hosting platform to set environment variables or secrets.
