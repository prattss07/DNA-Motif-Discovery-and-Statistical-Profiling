"""
DNA Motif Discovery & Clinical Profiling Dashboard:
Interactive Streamlit UI with Authentication (Sign Up, Sign In, Logout),
Live Z-Algorithm Genomic Scanner, Real-Time Repeat Profiling, Visual Charts,
and User Analytics History Tracking.
"""

import os
import sys
import pickle
import urllib.request
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# Add current and parent directory to sys.path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(BASE_DIR, "src")
DATA_DIR = os.path.join(BASE_DIR, "data")
MODELS_DIR = os.path.join(BASE_DIR, "models")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from auth import register_user, authenticate_user, record_user_analysis, get_user_analyses
from disease_registry import DISEASE_REGISTRY
from analyzer import read_fasta
from feature_extractor import extract_z_features
from train_model import FEATURE_COLUMNS, train_models

# Page Config
st.set_page_config(
    page_title="DNA Motif & Disease Profiling | Z-Algorithm",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling
st.markdown("""
<style>
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 8px;
        padding: 16px;
        border-left: 5px solid #0d6efd;
        margin-bottom: 12px;
    }
    .status-normal {
        color: #198754;
        font-weight: 700;
        font-size: 1.2rem;
    }
    .status-carrier {
        color: #ffc107;
        font-weight: 700;
        font-size: 1.2rem;
    }
    .status-pathogenic {
        color: #dc3545;
        font-weight: 700;
        font-size: 1.2rem;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        font-size: 16px;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)


# Session State Management
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


def render_auth_page():
    st.title("🧬 DNA Motif Discovery & Genetic Risk Profiling")
    st.subheader("High-Performance Genomic Sequence Analysis using the Z-Algorithm")

    col1, col2 = st.columns([1, 1], gap="large")

    with col1:
        st.markdown("### 🔐 Researcher / Clinician Sign In")
        with st.form("signin_form"):
            login_username = st.text_input("Username or Email")
            login_password = st.text_input("Password", type="password")
            submit_login = st.form_submit_button("Sign In", use_container_width=True)

            if submit_login:
                user = authenticate_user(login_username, login_password)
                if user:
                    st.session_state.authenticated = True
                    st.session_state.user = user
                    st.success(f"Welcome back, {user['full_name']}!")
                    st.rerun()
                else:
                    st.error("Invalid credentials. Please verify your username and password.")

    with col2:
        st.markdown("### 📝 New Account Registration")
        with st.form("signup_form"):
            reg_name = st.text_input("Full Name (e.g. Dr. John Doe)")
            reg_inst = st.text_input("Institution / Hospital (Optional)")
            reg_user = st.text_input("Choose Username")
            reg_email = st.text_input("Email Address")
            reg_pass = st.text_input("Choose Password", type="password")
            submit_signup = st.form_submit_button("Sign Up", use_container_width=True)

            if submit_signup:
                ok, msg = register_user(reg_user, reg_email, reg_pass, reg_name, reg_inst)
                if ok:
                    st.success(msg)
                else:
                    st.error(msg)


def render_dashboard():
    user = st.session_state.user

    # Sidebar Header & User Profile
    with st.sidebar:
        st.markdown(f"### 👤 Logged in as:")
        st.markdown(f"**{user['full_name']}**")
        if user['institution']:
            st.caption(f"🏛️ {user['institution']}")
        st.caption(f"📧 {user['email']}")

        if st.button("🚪 Log Out", use_container_width=True):
            st.session_state.authenticated = False
            st.session_state.user = None
            st.rerun()

        st.divider()
        st.markdown("### 🎯 Target Disease")
        selected_disease_name = st.selectbox(
            "Select Condition",
            options=list(DISEASE_REGISTRY.keys()),
            format_func=lambda x: x.replace("_", " ")
        )
        disease_info = DISEASE_REGISTRY[selected_disease_name]

        st.info(f"""
        **Target Gene:** `{disease_info['gene']}`  
        **Motif Pattern:** `{disease_info['motif']}`  
        **Chromosome:** `{disease_info['chromosome']}`  
        **Normal Limit:** $\le {disease_info['clinical_thresholds']['normal_max']}$ repeats  
        **Pathogenic:** $\ge {disease_info['clinical_thresholds']['pathogenic_min']}$ repeats
        """)

    # Main Navigation Tabs
    tab_scan, tab_user_analytics, tab_benchmark, tab_about = st.tabs([
        "🔬 Live Genomic Motif Scanner",
        "📊 User Analysis History & Analytics",
        "📁 Real NCBI Dataset Benchmark",
        "📖 DAA & Clinical Methodology"
    ])

    # TAB 1: LIVE GENOMIC SCANNER
    with tab_scan:
        st.header(f"Analyze DNA for {selected_disease_name.replace('_', ' ')}")
        st.markdown("Submit an authentic sequence using any of the 3 input methods below:")

        input_choice = st.radio(
            "Choose Sequence Source:",
            ["Choose from Downloaded NCBI GenBank Samples", "Fetch Live NCBI Accession via E-utilities", "Upload Local FASTA File", "Enter Raw DNA String"],
            horizontal=True
        )

        sequence = ""
        sample_name = ""

        if input_choice == "Choose from Downloaded NCBI GenBank Samples":
            fasta_files = [f for f in os.listdir(DATA_DIR) if f.endswith(".fasta")]
            # Filter files related to selected disease if possible
            matching_files = [f for f in fasta_files if disease_info['gene'] in f or disease_info['motif'] in f or selected_disease_name in f]
            display_files = matching_files if matching_files else fasta_files

            selected_file = st.selectbox("Select Specimen", display_files)
            if selected_file:
                sample_name = selected_file
                header, sequence = read_fasta(os.path.join(DATA_DIR, selected_file))
                st.caption(f"**FASTA Header:** `{header}`")

        elif input_choice == "Fetch Live NCBI Accession via E-utilities":
            col_acc, col_btn = st.columns([3, 1])
            with col_acc:
                acc_input = st.text_input("Enter NCBI Accession Number (e.g. KJ535072.1, NM_002111, NM_002024, U43748)", value=disease_info["ncbi_accessions"][0])
            with col_btn:
                st.write("")
                st.write("")
                fetch_btn = st.button("Fetch from NCBI", use_container_width=True)

            if fetch_btn and acc_input:
                with st.spinner(f"Querying NCBI Entrez API for {acc_input}..."):
                    try:
                        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc_input.strip()}&rettype=fasta&retmode=text"
                        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                        with urllib.request.urlopen(req) as resp:
                            fasta_content = resp.read().decode('utf-8')
                        save_path = os.path.join(DATA_DIR, f"{acc_input.strip()}.fasta")
                        with open(save_path, "w", encoding="utf-8") as f:
                            f.write(fasta_content)
                        header, sequence = read_fasta(save_path)
                        sample_name = f"NCBI Accession {acc_input.strip()}"
                        st.success(f"Downloaded {len(sequence):,} bp from NCBI GenBank!")
                    except Exception as e:
                        st.error(f"Error fetching from NCBI: {e}")

        elif input_choice == "Upload Local FASTA File":
            uploaded_file = st.file_uploader("Upload .fasta or .fa file", type=["fasta", "fa"])
            if uploaded_file is not None:
                content = uploaded_file.read().decode("utf-8")
                lines = content.splitlines()
                header = lines[0][1:] if lines[0].startswith(">") else "Uploaded Sample"
                sequence = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
                sample_name = uploaded_file.name

        else:
            raw_input = st.text_area("Paste Raw DNA Nucleotides (A, C, G, T)", height=120)
            if raw_input:
                sequence = "".join(raw_input.split()).upper()
                sample_name = "Manual Sequence Entry"

        # Action: Execute Z-Algorithm Analysis
        if sequence:
            st.divider()
            if st.button("🚀 Run Z-Algorithm & Predict Clinical Risk", type="primary", use_container_width=True):
                with st.spinner("Executing Linear-Time Z-Algorithm across DNA sequence..."):
                    # 1. Feature Extraction via Z-Algorithm
                    motif = disease_info["motif"]
                    gene = disease_info["gene"]
                    thresh = disease_info["clinical_thresholds"]

                    feats = extract_z_features(sequence, motif, gene)
                    repeats = feats["max_consecutive_repeats"]

                    # 2. ML Inference
                    clf = load_classifier()
                    df_feat = pd.DataFrame([feats])[FEATURE_COLUMNS]
                    pred_class = clf.predict(df_feat)[0]
                    pred_probs = clf.predict_proba(df_feat)[0]

                    label_names = ["Normal (Healthy)", "Carrier (Intermediate)", "Reduced Penetrance", "Pathogenic (Positive)"]
                    predicted_label = label_names[pred_class]
                    confidence = pred_probs[pred_class] * 100.0

                    # 3. Clinical Deterministic Status
                    if repeats >= thresh["pathogenic_min"]:
                        clinical_status = "PATHOGENIC / FULL PENETRANCE"
                        risk_level = "High"
                        status_css = "status-pathogenic"
                    elif repeats > thresh["normal_max"]:
                        clinical_status = "INTERMEDIATE / CARRIER RISK"
                        risk_level = "Moderate"
                        status_css = "status-carrier"
                    else:
                        clinical_status = "NORMAL (WITHIN REFERENCE LIMITS)"
                        risk_level = "Low"
                        status_css = "status-normal"

                    # 4. Save to User Analytics DB
                    record_user_analysis(user["id"], {
                        "sample_name": sample_name,
                        "disease": selected_disease_name,
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

                    st.markdown("### 📋 Clinical Diagnostic Summary")
                    c1, c2, c3, c4 = st.columns(4)
                    c1.metric("Sequence Length", f"{len(sequence):,} bp")
                    c2.metric(f"Total '{motif}' Hits", f"{feats['total_motif_occurrences']}")
                    c3.metric("Max Consecutive Repeats", f"{repeats} repeats")
                    c4.metric("Z-Score Enrichment", f"{feats['z_enrichment_score']} σ")

                    # Status Banner
                    st.markdown(f"""
                    <div style="padding:15px; border-radius:8px; background-color:#eef2f7; margin: 15px 0;">
                        <span style="font-size: 1.1rem; color: #495057;">Diagnosis:</span>
                        <div class="{status_css}">{clinical_status} [{risk_level} Risk]</div>
                        <div style="color: #6c757d; font-size: 0.95rem; margin-top: 5px;">
                            Machine Learning Classification: <b>{predicted_label}</b> (Confidence: <b>{confidence:.1f}%</b>)
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                    # Two Columns: DAA Telemetry vs Biological Profile
                    col_daa, col_chart = st.columns([1, 1], gap="medium")

                    with col_daa:
                        st.markdown("#### ⚡ DAA Algorithmic Telemetry (Z-Algorithm)")
                        st.write(f"- **Direct Character Comparisons:** `{feats['z_comparisons']:,}`")
                        st.write(f"- **Redundant Checks Bypassed via [L, R] Box:** `{feats['comparisons_saved']:,}`")
                        st.write(f"- **Algorithm Efficiency:** `{feats['efficiency_ratio'] * 100:.2f}%` checks saved vs naive scan")
                        st.write(f"- **Detected Repeat Locus:** Nucleotide coordinates `[{feats['locus_start']}:{feats['locus_end']}]`")

                        # Gauge Bar of comparisons saved
                        fig_gauge = go.Figure(go.Indicator(
                            mode="gauge+number",
                            value=feats['efficiency_ratio'] * 100,
                            title={'text': "Z-Algorithm Efficiency (% Checks Saved)"},
                            gauge={'axis': {'range': [0, 100]},
                                   'bar': {'color': "#0d6efd"},
                                   'steps': [
                                       {'range': [0, 50], 'color': "#e9ecef"},
                                       {'range': [50, 100], 'color': "#d1e7dd"}
                                   ]}
                        ))
                        fig_gauge.update_layout(height=240, margin=dict(l=20, r=20, t=30, b=20))
                        st.plotly_chart(fig_gauge, use_container_width=True)

                    with col_chart:
                        st.markdown("#### 🧠 ML Class Probability Distribution")
                        prob_df = pd.DataFrame({
                            "Risk Class": label_names,
                            "Probability": pred_probs * 100.0
                        })
                        fig_prob = px.bar(
                            prob_df,
                            x="Probability",
                            y="Risk Class",
                            orientation='h',
                            text="Probability",
                            color="Probability",
                            color_continuous_scale="Blues",
                            range_x=[0, 100]
                        )
                        fig_prob.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
                        fig_prob.update_layout(height=240, margin=dict(l=20, r=20, t=20, b=20), showlegend=False)
                        st.plotly_chart(fig_prob, use_container_width=True)

    # TAB 2: USER ANALYSIS HISTORY & ANALYTICS
    with tab_user_analytics:
        st.header(f"📊 Activity & Diagnostic History for {user['full_name']}")
        history = get_user_analyses(user["id"])

        if not history:
            st.info("No analyses recorded yet. Run a sequence in the 'Live Genomic Motif Scanner' tab to populate your history.")
        else:
            df_hist = pd.DataFrame(history)

            # Metric Counters
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Total Sequences Analyzed", len(df_hist))
            m2.metric("High-Risk Detections", len(df_hist[df_hist["risk_level"] == "High"]))
            total_checks_saved = df_hist["comparisons_saved"].sum()
            m3.metric("Total Algorithm Comparisons Saved", f"{total_checks_saved:,}")
            avg_repeats = df_hist["max_repeats"].mean()
            m4.metric("Avg Detected Repeat Length", f"{avg_repeats:.1f} repeats")

            st.divider()

            # Interactive Charts
            c_left, c_right = st.columns(2)

            with c_left:
                st.subheader("Disease Risk Breakdown")
                risk_counts = df_hist["risk_level"].value_counts().reset_index()
                risk_counts.columns = ["Risk Level", "Count"]
                fig_pie = px.pie(
                    risk_counts,
                    names="Risk Level",
                    values="Count",
                    color="Risk Level",
                    color_discrete_map={"Low": "#198754", "Moderate": "#ffc107", "High": "#dc3545"}
                )
                fig_pie.update_layout(height=300, margin=dict(l=10, r=10, t=10, b=10))
                st.plotly_chart(fig_pie, use_container_width=True)

            with c_right:
                st.subheader("Max Repeat Length by Sample")
                fig_bar = px.bar(
                    df_hist,
                    x="sample_name",
                    y="max_repeats",
                    color="risk_level",
                    color_discrete_map={"Low": "#198754", "Moderate": "#ffc107", "High": "#dc3545"},
                    labels={"sample_name": "Sample", "max_repeats": "Consecutive Repeats"}
                )
                fig_bar.update_layout(height=300, margin=dict(l=10, r=10, t=10, b=10))
                st.plotly_chart(fig_bar, use_container_width=True)

            st.subheader("Detailed Audit Log")
            st.dataframe(
                df_hist[[
                    "id", "sample_name", "disease", "target_motif",
                    "seq_length", "max_repeats", "prediction_label",
                    "risk_level", "confidence_pct", "created_at"
                ]],
                use_container_width=True
            )

    # TAB 3: REAL NCBI DATASET BENCHMARK
    with tab_benchmark:
        st.header("📁 Peer-Reviewed NCBI GenBank Dataset Benchmark")
        st.markdown("""
        All sequences below are 100% genuine biological specimens deposited in NCBI GenBank 
        by medical research laboratories.
        """)

        json_summary_path = os.path.join(BASE_DIR, "results", "real_ncbi_analysis_summary.json")
        if os.path.exists(json_summary_path):
            with open(json_summary_path, "r", encoding="utf-8") as jf:
                import json
                summary_data = json.load(jf)
            
            records_table = []
            for item in summary_data:
                feats = item["features"]
                records_table.append({
                    "Accession": item["accession"],
                    "Target Gene": item["gene"],
                    "Disease": item["disease"].replace("_", " "),
                    "Motif": item["target_motif"],
                    "Length (bp)": f"{item['sequence_length_bp']:,}",
                    "Max Repeats": feats["max_consecutive_repeats"],
                    "Z-Comparisons": f"{feats['z_comparisons']:,}",
                    "Saved Comparisons": f"{feats['comparisons_saved']:,}",
                    "Skip Efficiency": f"{feats['efficiency_ratio']*100:.1f}%",
                    "Clinical Status": item["clinical_status"]
                })

            st.dataframe(pd.DataFrame(records_table), use_container_width=True)
        else:
            st.info("Run `python src/run_real_pipeline.py` to generate the complete benchmark JSON.")

    # TAB 4: ABOUT & METHODOLOGY
    with tab_about:
        st.header("📖 Algorithmic & Bioinformatic Foundations")
        st.markdown(r"""
        ### 1. The Classical Z-Algorithm
        For a string $S$ of length $N$, the Z-array stores at each index $i > 0$ the length of the longest substring 
        starting at $S[i]$ that matches a prefix of $S$:
        $$Z[i] = \max \{ k \mid S[i \dots i+k-1] = S[0 \dots k-1] \}$$

        By maintaining an active window $[L, R]$ (the rightmost match segment), the algorithm achieves strictly 
        **linear time $O(|P| + |T|)$** by avoiding redundant character comparisons.

        ### 2. Clinical Triplet Repeat Genetics
        Trinucleotide repeat disorders are neurological conditions characterized by uncontrollable expansions 
        of 3-nucleotide motifs within specific human genes. The American College of Medical Genetics (ACMG) 
        defines standardized thresholds:
        * **Huntington's Disease (*HTT* / `CAG`):** Normal $\le 26$, Carrier $27\text{--}35$, Pathogenic $\ge 40$.
        * **Fragile X Syndrome (*FMR1* / `CGG`):** Normal $\le 44$, Pre-mutation $55\text{--}200$, Full Mutation $> 200$.
        * **Friedreich's Ataxia (*FXN* / `GAA`):** Normal $\le 33$, Carrier $34\text{--}65$, Pathogenic $\ge 100$.
        """)


def main():
    if not st.session_state.authenticated:
        render_auth_page()
    else:
        render_dashboard()


if __name__ == "__main__":
    main()
