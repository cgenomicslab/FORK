import streamlit as st

st.set_page_config(page_title="UniProt Lab Manager", layout="wide")


st.image("logo.png", width=300) 

st.title("UniProt Lab Manager")
st.markdown("---")
st.markdown("""
Welcome to the CGLab UniProt Reference Proteomes tool.

**Select a module from the sidebar** to get started.

### Available modules
- **UniProt Lab Manager** — sequence retrieval, HMM search, phylogenetic trees, presence/absence analysis
""")
