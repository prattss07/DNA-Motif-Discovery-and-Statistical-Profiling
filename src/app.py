"""
DNA Motif Discovery & Clinical Diagnostic Suite
Clean, simple user interface with standard structured clinical reporting
and an explicit Run button for pasted DNA sequences.
"""

import os
import sys
import pickle
import json
from datetime import datetime
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
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', sans-serif;
    }

    code, pre {
        font-family: 'JetBrains Mono', monospace !important;
    }

    /* Standard Medical Report Box */
    .clinical-report-card {
        border: 1px solid rgba(128, 128, 128, 0.25);
        border-radius: 10px;
        padding: 22px 26px;
        background-color: rgba(128, 128, 128, 0.04);
        margin: 20px 0;
    }

    .report-table {
        width: 100%;
        border-collapse: collapse;
        margin: 12px 0;
    }

    .report-table td {
        padding: 8px 12px;
        border-bottom: 1px solid rgba(128, 128, 128, 0.15);
        font-size: 0.95rem;
    }

    .report-table td.field-label {
        font-weight: 600;
        color: #64748b;
        width: 32%;
    }

    .report-table td.field-value {
        font-weight: 500;
    }

    /* Diagnosis Badges */
    .badge-normal {
        color: #10b981;
        font-weight: 700;
        font-size: 1.15rem;
    }
    .badge-carrier {
        color: #f59e0b;
        font-weight: 700;
        font-size: 1.15rem;
    }
    .badge-pathogenic {
        color: #ef4444;
        font-weight: 700;
        font-size: 1.15rem;
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
        should_run = False

        if input_type == "NCBI GenBank Clones":
            fasta_files = [f for f in os.listdir(DATA_DIR) if f.endswith(".fasta")]
            matched = [f for f in fasta_files if d_info['gene'] in f or selected_disease in f]
            options = matched if matched else fasta_files
            selected_file = st.selectbox("Select Specimen File", options)
            if selected_file:
                sample_name = selected_file
                header, sequence = read_fasta(os.path.join(DATA_DIR, selected_file))
                st.caption(f"Sequence ID: `{header[:85]}`")
                if st.button("Run Diagnostic Analysis", type="primary", use_container_width=True):
                    should_run = True

        elif input_type == "Upload FASTA File":
            up = st.file_uploader("Upload .fasta or .fa file", type=["fasta", "fa"])
            if up:
                lines = up.read().decode("utf-8").splitlines()
                header = lines[0][1:] if lines[0].startswith(">") else "Custom Upload"
                sequence = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
                sample_name = up.name
                if st.button("Run Diagnostic Analysis", type="primary", use_container_width=True):
                    should_run = True

        else:
            # Explicit Form for Pasting DNA String to provide an immediate Run button without Ctrl+Enter
            with st.form("paste_dna_form"):
                raw_seq = st.text_area("Paste DNA Sequence (A, C, G, T)", height=140, placeholder="ATGGCGACCCTGGAAAAGTCCTTCCAGCAGCAGCAGCAGCAGCAG...")
                run_pasted = st.form_submit_button("Run Diagnostic Analysis", type="primary", use_container_width=True)

                if run_pasted:
                    if raw_seq and raw_seq.strip():
                        sequence = "".join(raw_seq.split()).upper()
                        sample_name = "Pasted Nucleotides"
                        should_run = True
                    else:
                        st.warning("Please paste a DNA nucleotide string to analyze.")

        # Execute Diagnostic Evaluation
        if should_run and sequence:
            with st.spinner("Executing sequence verification..."):
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

                # 3. Clinical Status & Badge
                if repeats >= thresh["pathogenic_min"]:
                    clinical_status = "PATHOGENIC / FULL PENETRANCE"
                    risk_level = "High Risk"
                    badge_class = "badge-pathogenic"
                    recommendation = "Positive for pathological expansion. Referral to specialized neurogenetics clinic and supportive multidisciplinary care advised."
                elif repeats > thresh["normal_max"]:
                    clinical_status = "INTERMEDIATE / CARRIER PREMUTATION"
                    risk_level = "Moderate Risk"
                    badge_class = "badge-carrier"
                    recommendation = "Intermediate allele detected. Patient is asymptomatic; genetic counseling recommended for reproductive and transmission risk to offspring."
                else:
                    clinical_status = "NORMAL / WITHIN REFERENCE LIMITS"
                    risk_level = "Low Risk"
                    badge_class = "badge-normal"
                    recommendation = "Normal stable allele. No pathological expansion detected; no medical follow-up indicated."

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

                now_date = datetime.now().strftime("%B %d, %Y - %H:%M")

                # STANDARD STRUCTURED CLINICAL REPORT
                st.markdown(f"""
                <div class="clinical-report-card">
                    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid rgba(128,128,128,0.2); padding-bottom:8px; margin-bottom:14px;">
                        <span style="font-weight:700; font-size:1.1rem; text-transform:uppercase; letter-spacing:0.04em;">Molecular Genetics Pathology Report</span>
                        <span style="color:#64748b; font-size:0.85rem;">Date: {now_date}</span>
                    </div>

                    <table class="report-table">
                        <tr>
                            <td class="field-label">Specimen Identifier</td>
                            <td class="field-value"><code>{sample_name}</code> ({len(sequence):,} base pairs)</td>
                        </tr>
                        <tr>
                            <td class="field-label">Condition & Gene Panel</td>
                            <td class="field-value">{selected_disease.replace('_', ' ')} (<code>{gene}</code> Gene - Chromosome {d_info['chromosome']})</td>
                        </tr>
                        <tr>
                            <td class="field-label">Target Nucleotide Motif</td>
                            <td class="field-value"><code>5'-{motif}-3'</code> ({feats['total_motif_occurrences']} total matches identified)</td>
                        </tr>
                        <tr>
                            <td class="field-label">Identified Repeat Tract</td>
                            <td class="field-value"><b>{repeats} Consecutive Repeats</b> (Locus: nucleotide {feats['locus_start']} to {feats['locus_end']})</td>
                        </tr>
                        <tr>
                            <td class="field-label">Diagnostic Finding</td>
                            <td class="field-value"><span class="{badge_class}">{clinical_status}</span> &nbsp;({risk_level})</td>
                        </tr>
                        <tr>
                            <td class="field-label">Machine Learning Classification</td>
                            <td class="field-value">{predicted_label} (Confidence: <b>{confidence:.1f}%</b>)</td>
                        </tr>
                        <tr>
                            <td class="field-label">Clinical Action / Recommendation</td>
                            <td class="field-value">{recommendation}</td>
                        </tr>
                    </table>
                </div>
                """, unsafe_allow_html=True)

                # Summary Counters
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Sequence Length", f"{len(sequence):,} bp")
                c2.metric("Max Repeat Count", f"{repeats} repeats")
                c3.metric(f"Total '{motif}' Matches", f"{feats['total_motif_occurrences']}")
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
                color_discrete_map={"Low Risk": "#10b981", "Moderate Risk": "#f59e0b", "High Risk": "#ef4444",
                                    "Low": "#10b981", "Moderate": "#f59e0b", "High": "#ef4444"},
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
