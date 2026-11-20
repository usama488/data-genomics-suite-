"""Chest X-ray AI — Streamlit app.  Run: streamlit run app.py"""
import streamlit as st
import mx_page_learn as page_learn  # noqa: E402
import mx_page_study as page_study  # noqa: E402
import mx_page_worklist as page_worklist  # noqa: E402

PAGES = {"🩻 Study viewer": page_study.render, "📋 Triage worklist": page_worklist.render, "📚 Learn & model card": page_learn.render}


def main():
    st.sidebar.title("🩻 Chest X-ray AI")
    page = st.sidebar.radio("Navigate", list(PAGES), label_visibility="collapsed")
    s = st.session_state.get("mx_study")
    if s:
        st.sidebar.caption(f"Open study: {s['name']} · triage {s['res']['triage']}")
    st.sidebar.caption("⚠️ Research/education only — not a medical device.")
    PAGES[page]()


