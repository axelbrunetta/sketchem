import streamlit as st
import os
from pathlib import Path
import base64
from sketchem.utils.back_button import back_button
from sketchem.utils.environment import is_running_locally


def render_guide_page():
<<<<<<< HEAD
    # Add custom CSS for the page layout
    st.markdown(
        """
        <style>
        /* Hide the fullscreen button */
        .css-1lb4qv9 {display: none;}
        
        /* Adjust main container padding */
        .block-container {
            padding-top: 2rem;
            padding-bottom: 0rem;
            padding-left: 2rem;
            padding-right: 2rem;
        }
        
        /* PDF container styling */
        .pdf-container {
            display: flex;
            justify-content: center;
            width: 100%;
            margin-top: 1rem;
        }
        
        /* PDF iframe styling */
        .pdf-frame {
            width: 95%;
            height: 85vh;
            border: none;
            border-radius: 8px;
            box-shadow: 0 2px 6px rgba(0, 0, 0, 0.1);
        }

        /* Center the title */
        h1 {
            text-align: center !important;
            padding-bottom: 1rem;
        }
        </style>
    """,
        unsafe_allow_html=True,
    )

    st.title("A guide on how to use sketchem")

    # Construct path to the userguide.pdf in the data directory
    current_file = Path(__file__)
    src_dir = current_file.parent.parent
    pdf_path = os.path.join(src_dir, "data", "userguide.pdf")

=======
    # Import CSS from file
    css_path = os.path.join(os.path.dirname(__file__), "style", "guide_page_styling.css") if is_running_locally() else '/mount/src/sketchem/src/sketchem/pages/style/guide_page_styling.css'
    
    with open(css_path) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Add back button that redirects to home
    back_button(destination=None, label="Back to Home")
    
    st.markdown("<h2 style='text-align: center;'>A guide on how to use sketchem</h2>", unsafe_allow_html=True)
    
    

    # Get user's home directory and construct path to Downloads
    current_dir = os.path.dirname(__file__)
    pdf_path = os.path.join(current_dir, "..", "data", "userguide.pdf")
    
>>>>>>> origin/axel-s-branch-3
    # Check if the PDF file exists
    if os.path.exists(pdf_path):
        # Read and display the PDF file
        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()

            # Display PDF using PDF viewer with custom styling
            base64_pdf = base64.b64encode(pdf_bytes).decode("utf-8")
            pdf_display = f"""
                <div class="pdf-container">
                    <iframe 
                        src="data:application/pdf;base64,{base64_pdf}"
                        width="700" height="1000" 
                        class="pdf-frame"
                        type="application/pdf">
                    </iframe>
                </div>
            """
            st.markdown(pdf_display, unsafe_allow_html=True)
    else:
        st.error(
            f"User guide PDF file not found. Please make sure 'userguide.pdf' exists in the data directory at: {pdf_path}"
        )


if __name__ == "__main__":
    render_guide_page()