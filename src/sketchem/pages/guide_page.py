import streamlit as st
import base64
import os

# ---- TITLE AND SUBTITLE ----
st.markdown("<h1 style='text-align: center;'>🧪 Sketchem</h1>", unsafe_allow_html=True)
st.markdown("<h2 style='text-align: center;'>A guide on how to use sketchem</h2>", unsafe_allow_html=True)

# ---- DISPLAY PDF ----
try:
    current_dir = os.path.dirname(__file__)
    pdf_path = os.path.join(current_dir, "..", "static", "finalguide.pdf")

    with open(pdf_path, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode("utf-8")

    pdf_display = f'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
    st.markdown(pdf_display, unsafe_allow_html=True)

except FileNotFoundError:
    st.error("User guide PDF file not found.")