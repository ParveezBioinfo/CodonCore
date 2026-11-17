import streamlit as st
import plotly.graph_objects as go
from analyzer import get_basic_stats, calculate_sliding_window_gc, find_orfs

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="Codoncore", layout="wide")

# --- LEFT DARK SIDEBAR ---
st.sidebar.title("🧬 Codoncore")
st.sidebar.markdown("Bioinformatics Web Dashboard")

uploaded_file = st.sidebar.file_uploader("Upload FASTA File", type=["fasta", "fa"])
manual_sequence = st.sidebar.text_area("Or Paste Manual Sequence Input")

min_orf_length = st.sidebar.slider("Minimum ORF Length (Amino Acids)", min_value=10, max_value=300, value=50)

analyze_button = st.sidebar.button("Analyze Sequence", type="primary")

# --- RIGHT LIGHT MAIN DASHBOARD ---
st.title("Sequence Analysis Results")

if analyze_button:
    sequence_data = ""
    
    # Process the uploaded file or text input
    if uploaded_file is not None:
        raw_text = uploaded_file.getvalue().decode("utf-8").split("\n")
        sequence_data = "".join([line.strip() for line in raw_text if not line.startswith(">")])
    elif manual_sequence:
        sequence_data = "".join(manual_sequence.split())
        
    if not sequence_data:
        st.sidebar.error("Please upload a file or paste a sequence first!")
    else:
        # 1. Pull data from your analyzer.py
        seq_length, gc_content = get_basic_stats(sequence_data)
        orf_df = find_orfs(sequence_data, min_orf_length)
        
        # 2. Top Summary Cards
        col1, col2, col3 = st.columns(3)
        col1.metric("Sequence Length (bp)", f"{seq_length:,}")
        col2.metric("GC Content (%)", f"{gc_content}%")
        col3.metric("Detected ORFs", len(orf_df))
        
        st.markdown("---")
        
        # 3. Interactive Plotly Chart
        st.subheader("Sliding-Window GC Content (%)")
        positions, gc_values = calculate_sliding_window_gc(sequence_data)
        
        fig = go.Figure(data=go.Scatter(x=positions, y=gc_values, mode='lines', name='GC%'))
        fig.update_layout(xaxis_title="Position (bp)", yaxis_title="GC Content (%)", margin=dict(l=0, r=0, t=30, b=0))
        st.plotly_chart(fig, use_container_width=True)
        
        st.markdown("---")
        
        # 4. Clean Data Table
        st.subheader("6-Frame Open Reading Frames (ORFs)")
        if len(orf_df) > 0:
            st.dataframe(orf_df, use_container_width=True, hide_index=True)
        else:
            st.info("No ORFs found meeting the minimum length criteria.")