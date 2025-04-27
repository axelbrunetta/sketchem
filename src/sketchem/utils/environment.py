
import platform

# Behavior for streamlit is different in cloud and locally run -> especially for what we will do with the databse, hence the need for this

def is_running_locally():
    """Check if the app is running locally or in Streamlit cloud -> this works since Streamlit Cloud runs on Linux without a processor name"""
    return platform.processor() != ''
