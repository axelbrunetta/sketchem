import streamlit as st
import os
from pathlib import Path
import base64

def render_guide_page():
    # Add custom CSS for the page layout
    st.markdown("""
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
    """, unsafe_allow_html=True)
    
    st.title("A guide on how to use SketChem")
    
    # Get user's home directory and construct path to Downloads
    home = str(Path.home())
    pdf_path = os.path.join(home, "Downloads", "userguide_final.pdf")
    
    # Check if the PDF file exists
    if os.path.exists(pdf_path):
        # Read and display the PDF file
        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()
            
            # Display PDF using PDF viewer with custom styling
            base64_pdf = base64.b64encode(pdf_bytes).decode('utf-8')
            pdf_display = f'''
                <div class="pdf-container">
                    <iframe 
                        src="data:application/pdf;base64,{base64_pdf}" 
                        class="pdf-frame"
                        type="application/pdf">
                    </iframe>
                </div>
            '''
            st.markdown(pdf_display, unsafe_allow_html=True)
    else:
        st.error(f"User guide PDF file not found. Please make sure 'userguide_final.pdf' exists in your Downloads folder at: {pdf_path}")

if __name__ == "__main__":
    render_guide_page()
