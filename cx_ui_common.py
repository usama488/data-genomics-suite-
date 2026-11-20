import streamlit as st

import cx_workflow as W
from cx_learn import STEP_INFO

ICON = {"warn": "⚠️", "good": "✅", "info": "💡"}
STATUS_ICON = {"PASS": "🟢", "WARN": "🟠", "FAIL": "🔴"}


def init():
    st.session_state.setdefault("cx_case", None)
    st.session_state.setdefault("learn", True)


def case() -> W.Case | None:
    return st.session_state.get("cx_case")


def need_case(*steps: str) -> W.Case:
    c = case()
    if c is None or not c.fastq:
        st.info("👈 Start on **🧾 Case & FASTQ** (or run the one-click demo on 🏠 Home).")
        st.stop()
    missing = [s for s in steps if s not in c.done]
    if missing:
        st.info(f"Run these steps first (▶️ Run pipeline page): **{', '.join(missing)}**")
        st.stop()
    return c


def tip(*terms):
    if not st.session_state.get("learn"):
        return
    for t in terms:
        if t in STEP_INFO:
            exp, ex = STEP_INFO[t]
            st.caption(f"🎓 **{t}** — {exp}" + (f" *{ex}*" if ex else ""))


def findings(fs):
    for f in fs:
        {"warn": st.warning, "good": st.success, "info": st.info}[f["level"]](f["text"], icon=ICON[f["level"]])


def show(fig, key=None, **kw):
    if fig is not None:
        return st.plotly_chart(fig, key=key, **kw)
