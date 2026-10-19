"""
DNA Motif Discovery & Clinical Diagnostic Suite
Clean, modern, emoji-free user interface focused on seamless diagnostic prediction.
"""

import os
import sys
import pickle
import urllib.request
import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# Setup Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(BASE_DIR, "src")
DATA_DIR = os.path.join(BASE_DIR, "data")
MODELS_DIR = os.path.join(BASE_DIR, "models")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from auth import register_user, authenticate_user, record_user_analysis, get_user_analyses
from disease_registry import DISEASE_REGISTRY
from analyzer import read_fasta
from feature_extractor import extract_z_features
from train_model import FEATURE_COLUMNS, train_models

# Page Configurations
st.set_page_config(
    page_title="GenomiX | Clinical Motif Platform",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Clean, Modern, Native Dark/Light Adaptive CSS (No hardcoded whites that cause blank blocks)
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', sans-serif;
    }

    code, pre {
        font-family: 'JetBrains Mono', monospace !important;
    }

    /* Clean Card Container */
    .clean-card {
        border: 1px solid rgba(128, 128, 128, 0.2);
        border-radius: 12px;
        padding: 24px 28px;
        margin-bottom: 20px;
    }

    /* Diagnosis Status Banners */
    .status-box-normal {
        border-left: 5px solid #10b981;
        background-color: rgba(16, 185, 129, 0.08);
        border-radius: 8px;
        padding: 16px 20px;
        margin: 16px 0;
    }
    .status-box-carrier {
        border-left: 5px solid #f59e0b;
        background-color: rgba(245, 158, 11, 0.08);
        border-radius: 8px;
        padding: 16px 20px;
        margin: 16px 0;
    }
    .status-box-pathogenic {
        border-left: 5px solid #ef4444;
        background-color: rgba(239, 68, 68, 0.08);
        border-radius: 8px;
        padding: 16px 20px;
        margin: 16px 0;
    }

    /* Simple Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
        border-bottom: 1px solid rgba(128, 128, 128, 0.2);
    }
    .stTabs [data-baseweb="tab"] {
        padding: 10px 4px;
        font-size: 0.95rem;
        font-weight: 500;
    }
</style>
""", unsafe_allow_html=True)

# Session State
if "authenticated" not in st.session_state:
    st.session_state.authenticated = False
if "user" not in st.session_state:
    st.session_state.user = None


def load_classifier():
    model_path = os.path.join(MODELS_DIR, "z_motif_rf_model.pkl")
    if not os.path.exists(model_path):
        train_models()
    with open(model_path, "rb") as f:
        return pickle.load(f)


def render_auth():
    st.title("GenomiX | Clinical Motif Platform")
    st.caption("Genomic Sequence Verification and Multi-Disease Repeat Expansion Diagnostic System")

    auth_tab1, auth_tab2 = st.tabs(["Sign In", "Create Account"])

    with auth_tab1:
        col1, col2, col3 = st.columns([1, 1.8, 1])
        with col2:
            st.markdown("### Clinical Sign In")
            st.caption("Access your patient sequence diagnostic workspace.")
            with st.form("signin_form"):
                login_id = st.text_input("Username or Email", placeholder="e.g. jsmith266")
                login_pass = st.text_input("Password", type="password", placeholder="Enter password")
                submit_login = st.form_submit_button("Sign In", type="primary", use_container_width=True)

                if submit_login:
                    if not login_id or not login_pass:
                        st.warning("Please enter your username and password.")
                    else:
                        user = authenticate_user(login_id, login_pass)
                        if user:
                            st.session_state.authenticated = True
                            st.session_state.user = user
                            st.rerun()
                        else:
                            st.error("Invalid credentials. Please verify your details.")

    with auth_tab2:
        col1, col2, col3 = st.columns([1, 1.8, 1])
        with col2:
            st.markdown("### Create New Account")
            st.caption("Register for genomic sequence analysis and audit tracking.")
            with st.form("signup_form"):
                full_name = st.text_input("Full Name & Title", placeholder="e.g. Dr. Jane Smith")
                institution = st.text_input("Institution (Optional)", placeholder="e.g. Mayo Clinic")
                c1, c2 = st.columns(2)
                with c1:
                    new_user = st.text_input("Username", placeholder="e.g. jsmith266")
                with c2:
                    new_email = st.text_input("Email", placeholder="e.g. jsmith@hospital.org")
                new_pass = st.text_input("Password", type="password", placeholder="Minimum 6 characters")

                submit_signup = st.form_submit_button("Create Account", use_container_width=True)

                if submit_signup:
                    if not new_user or not new_email or not new_pass or not full_name:
                        st.warning("Please complete all required fields.")
                    else:
                        ok, msg = register_user(new_user, new_email, new_pass, full_name, institution)
                        if ok:
                            st.success(msg)
                        else:
                            st.error(msg)


def render_dashboard():
    user = st.session_state.user

    # Sidebar: User profile and Disease Selection
    with st.sidebar:
        st.subheader("Current User")
        st.write(f"**{user['full_name']}**")
        if user['institution']:
            st.caption(f"{user['institution']}")
        st.caption(f"{user['email']}")

        if st.button("Log Out", use_container_width=True):
            st.session_state.authenticated = False
            st.session_state.user = None
            st.rerun()

        st.divider()
        st.subheader("Target Disease")
        selected_disease = st.selectbox(
            "Select Condition Panel",
            options=list(DISEASE_REGISTRY.keys()),
            format_func=lambda x: x.replace("_", " ")
        )
        d_info = DISEASE_REGISTRY[selected_disease]

        st.markdown(f"""
        - **Target Gene:** `{d_info['gene']}`
        - **Motif Pattern:** `{d_info['motif']}`
        - **Chromosome:** `{d_info['chromosome']}`
        - **Normal Range:** $\le {d_info['clinical_thresholds']['normal_max']}$ repeats
        - **Pathogenic:** $\ge {d_info['clinical_thresholds']['pathogenic_min']}$ repeats
        """)

    # Main Area
    st.title(f"{selected_disease.replace('_', ' ')} Diagnostics")
    st.caption(f"Screening for repeat expansions in the {d_info['gene']} gene ({d_info['motif']} motif)")

    tab_scan, tab_history, tab_benchmark = st.tabs([
        "Sequence Scanner",
        "Analysis History",
        "NCBI Reference Cohort"
    ])

    # TAB 1: SCANNER
    with tab_scan:
        input_type = st.radio(
            "Sequence Input Method:",
            ["NCBI GenBank Clones", "Query Live Accession (NCBI)", "Upload FASTA File", "Paste DNA String"],
            horizontal=True
        )

        sequence = ""
        sample_name = ""

        if input_type == "NCBI GenBank Clones":
            fasta_files = [f for f in os.listdir(DATA_DIR) if f.endswith(".fasta")]
            matched = [f for f in fasta_files if d_info['gene'] in f or selected_disease in f]
            options = matched if matched else fasta_files
            selected_file = st.selectbox("Select Specimen File", options)
            if selected_file:
                sample_name = selected_file
                header, sequence = read_fasta(os.path.join(DATA_DIR, selected_file))
                st.caption(f"Sequence ID: `{header[:85]}`")

        elif input_type == "Query Live Accession (NCBI)":
            col_acc, col_btn = st.columns([3, 1])
            with col_acc:
                acc = st.text_input("Enter NCBI Accession Number", value=d_info["ncbi_accessions"][0])
            with col_btn:
                st.write("")
                st.write("")
                fetch = st.button("Fetch Sequence", use_container_width=True)

            if fetch and acc:
                with st.spinner("Connecting to NCBI GenBank..."):
                    try:
                        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc.strip()}&rettype=fasta&retmode=text"
                        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                        with urllib.request.urlopen(req) as resp:
                            fasta_raw = resp.read().decode('utf-8')
                        save_p = os.path.join(DATA_DIR, f"{acc.strip()}.fasta")
                        with open(save_p, "w", encoding="utf-8") as f:
                            f.write(fasta_raw)
                        header, sequence = read_fasta(save_p)
                        sample_name = f"NCBI Accession {acc.strip()}"
                        st.success(f"Downloaded {len(sequence):,} base pairs from NCBI.")
                    except Exception as e:
                        st.error(f"Error fetching from NCBI: {e}")

        elif input_type == "Upload FASTA File":
            up = st.file_uploader("Upload .fasta or .fa file", type=["fasta", "fa"])
            if up:
                lines = up.read().decode("utf-8").splitlines()
                header = lines[0][1:] if lines[0].startswith(">") else "Custom Upload"
                sequence = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
                sample_name = up.name

        else:
            raw_seq = st.text_area("Paste DNA Sequence (A, C, G, T)", height=100)
            if raw_seq:
                sequence = "".join(raw_seq.split()).upper()
                sample_name = "Pasted Nucleotides"

        # Action: Analyze Sequence
        if sequence:
            st.write(f"Loaded Sequence: **{sample_name}** ({len(sequence):,} base pairs)")
            if st.button("Run Diagnostic Analysis", type="primary", use_container_width=True):
                with st.spinner("Analyzing sequence..."):
                    motif = d_info["motif"]
                    gene = d_info["gene"]
                    thresh = d_info["clinical_thresholds"]

                    # 1. Feature Extraction via Z-Algorithm
                    feats = extract_z_features(sequence, motif, gene)
                    repeats = feats["max_consecutive_repeats"]

                    # 2. Machine Learning Classifier
                    clf = load_classifier()
                    df_feat = pd.DataFrame([feats])[FEATURE_COLUMNS]
                    pred_class = clf.predict(df_feat)[0]
                    pred_probs = clf.predict_proba(df_feat)[0]

                    labels = ["Normal (Healthy)", "Carrier (Intermediate)", "Reduced Penetrance", "Pathogenic (Full Penetrance)"]
                    predicted_label = labels[pred_class]
                    confidence = pred_probs[pred_class] * 100.0

                    # 3. Clinical Status
                    if repeats >= thresh["pathogenic_min"]:
                        clinical_status = "PATHOGENIC / FULL PENETRANCE"
                        risk_level = "High"
                        box_class = "status-box-pathogenic"
                    elif repeats > thresh["normal_max"]:
                        clinical_status = "INTERMEDIATE / CARRIER"
                        risk_level = "Moderate"
                        box_class = "status-box-carrier"
                    else:
                        clinical_status = "NORMAL (WITHIN REFERENCE LIMITS)"
                        risk_level = "Low"
                        box_class = "status-box-normal"

                    # 4. Save to DB
                    record_user_analysis(user["id"], {
                        "sample_name": sample_name,
                        "disease": selected_disease,
                        "target_motif": motif,
                        "seq_length": len(sequence),
                        "max_repeats": repeats,
                        "total_matches": feats["total_motif_occurrences"],
                        "z_comparisons": feats["z_comparisons"],
                        "comparisons_saved": feats["comparisons_saved"],
                        "prediction_label": predicted_label,
                        "risk_level": risk_level,
                        "confidence_pct": round(confidence, 2)
                    })

                    # Status Banner
                    st.markdown(f"""
                    <div class="{box_class}">
                        <div style="font-size: 0.8rem; font-weight: 600; text-transform: uppercase;">Diagnostic Result</div>
                        <div style="font-size: 1.4rem; font-weight: 700; margin: 4px 0;">{clinical_status} [{risk_level.upper()} RISK]</div>
                        <div>Found <b>{repeats} consecutive '{motif}' repeats</b> (Coordinates: {feats['locus_start']}:{feats['locus_end']}). Predicted: <b>{predicted_label}</b> ({confidence:.1f}% confidence).</div>
                    </div>
                    """, unsafe_allow_html=True)

                    # Simple Metrics
                    c1, c2, c3, c4 = st.columns(4)
                    c1.metric("Sequence Length", f"{len(sequence):,} bp")
                    c2.metric("Max Repeat Length", f"{repeats} repeats")
                    c3.metric(f"Total '{motif}' Hits", f"{feats['total_motif_occurrences']}")
                    c4.metric("Model Confidence", f"{confidence:.1f}%")

                    # Class Probabilities Bar Chart
                    st.subheader("Prediction Confidence")
                    df_p = pd.DataFrame({"Class": labels, "Probability": pred_probs * 100.0})
                    fig_p = px.bar(
                        df_p, x="Probability", y="Class", orientation='h',
                        text="Probability", range_x=[0, 100]
                    )
                    fig_p.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
                    fig_p.update_layout(height=240, margin=dict(l=10, r=10, t=10, b=10))
                    st.plotly_chart(fig_p, use_container_width=True)

    # TAB 2: HISTORY
    with tab_history:
        st.subheader("Your Analysis History")
        history = get_user_analyses(user["id"])

        if not history:
            st.info("No recorded analyses yet. Run a sequence in the scanner tab to populate this history.")
        else:
            df_h = pd.DataFrame(history)

            m1, m2, m3 = st.columns(3)
            m1.metric("Total Analyses", len(df_h))
            m2.metric("High-Risk Detections", len(df_h[df_h["risk_level"] == "High"]))
            m3.metric("Average Repeat Count", f"{df_h['max_repeats'].mean():.1f} repeats")

            st.divider()

            # Repeat Count Chart
            fig_hist = px.bar(
                df_h, x="sample_name", y="max_repeats", color="risk_level",
                color_discrete_map={"Low": "#10b981", "Moderate": "#f59e0b", "High": "#ef4444"},
                labels={"sample_name": "Sample", "max_repeats": "Consecutive Repeats"}
            )
            fig_hist.update_layout(height=280, margin=dict(l=10, r=10, t=10, b=10))
            st.plotly_chart(fig_hist, use_container_width=True)

            st.subheader("Recent Runs")
            st.dataframe(
                df_h[[
                    "id", "sample_name", "disease", "target_motif",
                    "seq_length", "max_repeats", "prediction_label",
                    "risk_level", "confidence_pct", "created_at"
                ]],
                use_container_width=True
            )

    # TAB 3: BENCHMARK
    with tab_benchmark:
        st.subheader("NCBI GenBank Reference Cohort")
        st.caption("Peer-reviewed biological sequences verified from NCBI GenBank.")

        summary_file = os.path.join(RESULTS_DIR, "real_ncbi_analysis_summary.json")
        if os.path.exists(summary_file):
            with open(summary_file, "r", encoding="utf-8") as jf:
                s_data = json.load(jf)

            rows = []
            for item in s_data:
                f = item["features"]
                rows.append({
                    "Accession": item["accession"],
                    "Gene": item["gene"],
                    "Condition": item["disease"].replace("_", " "),
                    "Motif": item["target_motif"],
                    "Length (bp)": f"{item['sequence_length_bp']:,}",
                    "Max Repeats": f["max_consecutive_repeats"],
                    "Diagnosis": item["clinical_status"]
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True)


def main():
    if not st.session_state.authenticated:
        render_auth()
    else:
        render_dashboard()


if __name__ == "__main__":
    main()
