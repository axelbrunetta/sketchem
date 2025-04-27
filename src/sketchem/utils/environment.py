import os

# Behavior for streamlit is different in cloud and locally run -> especially for what we will do with the databse, hence the need for this

def is_running_locally():
    """Check if the app is running locally or in Streamlit cloud"""
    return not os.environ.get("STREAMLIT_DEPLOYMENT", False)
