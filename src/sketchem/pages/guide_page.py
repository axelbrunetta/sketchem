import streamlit as st
import os
from pathlib import Path
import base64
from sketchem.utils.back_button import back_button
from sketchem.utils.environment import is_running_locally


def render_guide_page():
    #import css
    css_path = os.path.join(os.path.dirname(__file__), "style", "guide_page_styling.css") if is_running_locally() else '/mount/src/sketchem/src/sketchem/pages/style/guide_page_styling.css'
    
    with open(css_path) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    #add back button that leads to home
    back_button(destination=None, label="Back to Home")
    
    st.markdown("<h2 style='text-align: center;'>A guide on how to use sketchem</h2>", unsafe_allow_html=True)
    
    

    #get user's home directory and construct path to downloads
    current_dir = os.path.dirname(__file__)
    #get path to pdf file
    pdf_path = os.path.join(current_dir, "..", "..", "..", "assets", "userguide.pdf")
    
    #check if pdf file exists
    if os.path.exists(pdf_path):
        #read and display pdf file
        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()

            #display pdf using pdf viewer 
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
    #error if pdf file not found
    else:
        st.error(
            f"User guide PDF file not found. Please make sure 'userguide.pdf' exists in the data directory at: {pdf_path}"
        )


if __name__ == "__main__":
    render_guide_page()