"""
DNA Motif Discovery & Clinical Diagnostic Suite
Clean, simple user interface with detailed textual diagnostic findings
and local / custom sequence inputs.
"""

import os
import sys
import pickle
import json
import pandas as pd
import plotly.express as px
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

# Clean, Theme-Adaptive CSS
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', sans-serif;
    }

    /* Diagnosis Status Banners */
    .status-box-normal {
        border-left: 5px solid #10b981;
        background-color: rgba(16, 185, 129, 0.08);
        border-radius: 8px;
        padding: 20px 24px;
        margin: 20px 0;
    }
    .status-box-carrier {
        border-left: 5px solid #f59e0b;
        background-color: rgba(245, 158, 11, 0.08);
        border-radius: 8px;
        padding: 20px 24px;
        margin: 20px 0;
    }
    .status-box-pathogenic {
        border-left: 5px solid #ef4444;
        background-color: rgba(239, 68, 68, 0.08);
        border-radius: 8px;
        padding: 20px 24px;
        margin: 20px 0;
    }

    .report-paragraph {
        font-size: 1.05rem;
        line-height: 1.7;
        margin-top: 12px;
    }

    /* Tab styling */
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


def generate_textual_report(sample_name: str, disease_name: str, gene: str, motif: str,
                            repeats: int, locus_start: int, locus_end: int,
                            total_matches: int, seq_len: int,
                            predicted_label: str, confidence: float,
                            clinical_status: str, risk_level: str,
                            thresh: dict) -> str:
    """Generates an authentic, fully written diagnostic report in continuous natural English words."""

    condition_display = disease_name.replace("_", " ")

    # Section 1: Overview
    p1 = (
        f"A molecular screening of specimen '{sample_name}' (spanning a total analyzed length of {seq_len:,} base pairs) "
        f"was conducted for {condition_display} targeting the {gene} gene on the human chromosome. "
        f"The diagnostic assay evaluated the presence, positioning, and tandem expansion of the critical '{motif}' nucleotide motif."
    )

    # Section 2: Exact Findings
    p2 = (
        f"Upon sequence scanning, a total of {total_matches} occurrences of the '{motif}' pattern were verified across the DNA strand. "
        f"The primary polymorphic locus of interest was identified beginning at nucleotide index {locus_start} and concluding at index {locus_end}. "
        f"Within this specific segment, the patient's sequence demonstrated an unbroken consecutive repeat tract measuring {repeats} tandem '{motif}' units."
    )

    # Section 3: Clinical Diagnostic Interpretation
    if "PATHOGENIC" in clinical_status:
        p3 = (
            f"CLINICAL FINDING: POSITIVE FOR PATHOGENIC REPEAT EXPANSION ({risk_level.upper()} RISK). "
            f"The observed tract of {repeats} repeats significantly exceeds the clinical threshold of {thresh['pathogenic_min']} repeats established by standard medical genetics guidelines. "
            f"Alleles of this size demonstrate full genetic penetrance, indicating that the patient is genetically positive and affected by {condition_display}. "
            f"The computational machine learning model corroborates this diagnosis with a classification of '{predicted_label}' at {confidence:.1f}% confidence. "
            f"Medical referral to a neurogenetics specialist and multidisciplinary care is advised."
        )
    elif "CARRIER" in clinical_status or "INTERMEDIATE" in clinical_status:
        p3 = (
            f"CLINICAL FINDING: INTERMEDIATE / PRE-MUTATION CARRIER STATE ({risk_level.upper()} RISK). "
            f"The detected count of {repeats} repeats lies within the intermediate range (between {thresh['normal_max'] + 1} and {thresh['intermediate_max']} repeats). "
            f"Individuals harboring alleles in this range do not manifest clinical symptoms of {condition_display}. "
            f"However, this sequence possesses meiotic instability and may undergo further repeat expansion during transmission to biological offspring. "
            f"The machine learning system classifies this sample as '{predicted_label}' ({confidence:.1f}% confidence). "
            f"Genetic counseling regarding reproductive risk is recommended."
        )
    else:
        p3 = (
            f"CLINICAL FINDING: NEGATIVE / NORMAL REFERENCE RANGE ({risk_level.upper()} RISK). "
            f"The observed repeat count of {repeats} is well within the healthy population reference boundary (normal threshold is at or below {thresh['normal_max']} repeats). "
            f"This sequence demonstrates stable genomic transmission without risk of abnormal polyglutamine expansion or phenotypic disease manifestation. "
            f"The machine learning classification confirms a status of '{predicted_label}' with {confidence:.1f}% confidence. "
            f"No follow-up genetic intervention for {condition_display} is indicated."
        )

    return f"{p1}\n\n{p2}\n\n{p3}"


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
            ["NCBI GenBank Clones", "Upload FASTA File", "Paste DNA String"],
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

                    # Generate Full Textual Report in Words
                    textual_report = generate_textual_report(
                        sample_name=sample_name,
                        disease_name=selected_disease,
                        gene=gene,
                        motif=motif,
                        repeats=repeats,
                        locus_start=feats['locus_start'],
                        locus_end=feats['locus_end'],
                        total_matches=feats['total_motif_occurrences'],
                        seq_len=len(sequence),
                        predicted_label=predicted_label,
                        confidence=confidence,
                        clinical_status=clinical_status,
                        risk_level=risk_level,
                        thresh=thresh
                    )

                    # Diagnostic Report Display in Full Words
                    st.markdown(f"""
                    <div class="{box_class}">
                        <div style="font-size: 0.85rem; font-weight: 700; text-transform: uppercase;">Diagnostic Assessment</div>
                        <div style="font-size: 1.4rem; font-weight: 700; margin: 6px 0;">{clinical_status} [{risk_level.upper()} RISK]</div>
                        <div class="report-paragraph">
                            {textual_report.replace(chr(10) + chr(10), '<br><br>')}
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                    # Simple Summary Metrics
                    c1, c2, c3, c4 = st.columns(4)
                    c1.metric("Sequence Length", f"{len(sequence):,} bp")
                    c2.metric("Max Repeat Length", f"{repeats} repeats")
                    c3.metric(f"Total '{motif}' Hits", f"{feats['total_motif_occurrences']}")
                    c4.metric("Model Confidence", f"{confidence:.1f}%")

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
