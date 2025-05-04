import streamlit as st

def fancy_gradient_button(label, key=None, color1="#FF4B4B", color2="#FFB47B", text_color="white"):
    """
    Creates a Streamlit button with a gradient background and custom text color
    using the robust Wrapper Div method.

    Args:
        label (str): The text label for the button.
        key (str, optional): An optional key for the Streamlit button widget.
                             Important for distinguishing button states if multiple
                             buttons are used. Defaults to None.
        color1 (str, optional): The starting color of the gradient (hex code).
                                Defaults to "#FF4B4B" (Streamlit red).
        color2 (str, optional): The ending color of the gradient (hex code).
                                Defaults to "#FFB47B" (Orangeish).
        text_color (str, optional): The color of the button text (CSS color name or hex code).
                                    Defaults to "white".

    Returns:
        bool: True if the button was clicked in the last run, False otherwise.
              Behaves exactly like the standard st.button().
    """
    button_id = f"gradient_button_{key or label.replace(' ', '_').lower()}"

    button_css = f"""
    <style>
    .fancy-gradient-button-container.{button_id} div[data-testid="stButton"] > button {{
        /* Gradient Background */
        background-image: linear-gradient(to right, {color1} 0%, {color2} 51%, {color1} 100%) !important;
        
        /* Text Styling */
        color: {text_color} !important;
        font-size: 16px !important;        /* Adjust font size */
        font-weight: bold !important;      /* Make text bold */
        text-align: center !important;     /* Center text */
        text-decoration: none !important;  /* Remove underline */

        /* Button Shape and Size */
        padding: 12px 25px !important;     /* Top/bottom and left/right padding */
        border-radius: 10px !important;    /* Rounded corners */
        border: none !important;           /* Remove default border */
        display: inline-block !important;  /* Allow padding and width */

        /* Interaction */
        cursor: pointer !important;        /* Hand cursor on hover */
        transition: 0.5s !important;       /* Smooth transition for hover effects */
        background-size: 200% auto !important; /* Crucial for the gradient shift hover effect */

        /* Optional: Add a shadow for depth */
        box-shadow: 0 4px 15px 0 rgba(116, 79, 168, 0.75) !important; /* Adjust color and spread as needed */
    }}

    /* Hover Effect: Shift the background gradient */
    .fancy-gradient-button-container.{button_id} div[data-testid="stButton"] > button:hover {{
        background-position: right center !important; /* Move gradient position */
        color: {text_color} !important;             /* Ensure text color stays */
        text-decoration: none !important;         /* Ensure no underline on hover */
        /* Optional: Enhance shadow on hover */
        box-shadow: 0 6px 20px 0 rgba(116, 79, 168, 0.85) !important;
    }}

    /* Active/Click Effect: Slight press down */
    .fancy-gradient-button-container.{button_id} div[data-testid="stButton"] > button:active {{
        transform: translateY(1px) !important; /* Move button down slightly */
        /* Optional: Reduce shadow on press */
        box-shadow: 0 2px 10px 0 rgba(116, 79, 168, 0.65) !important;
    }}
    </style>
    """
    st.markdown(button_css, unsafe_allow_html=True)
    st.markdown(f'<div class="fancy-gradient-button-container {button_id}">', unsafe_allow_html=True)

    clicked = st.button(label, key=key)

    st.markdown('</div>', unsafe_allow_html=True)
    return clicked
