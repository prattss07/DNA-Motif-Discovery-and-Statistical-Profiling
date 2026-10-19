"""
DNA Motif Discovery & Clinical Profiling Enterprise Suite
==========================================================
High-Performance Genomic Sequence Verification using the Linear-Time Z-Algorithm O(N+M)
Equipped with Professional Authentication, Real-Time Telemetry,
Comparative Visualizations, Audit Trails, and Clinical Diagnostics.
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

# Setup Pathing
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
    page_title="GenomiX | Clinical Motif Intelligence",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Enterprise Modern UI Styling
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap');

    html, body, [class*="css"] {
        font-family: 'Plus Jakarta Sans', sans-serif;
    }
    
    code, pre {
        font-family: 'JetBrains Mono', monospace !important;
    }

    /* Gradient Headers & Modern Cards */
    .hero-banner {
        background: linear-gradient(135deg, #0f172a 0%, #1e293b 50%, #0f766e 100%);
        color: white;
        padding: 38px 45px;
        border-radius: 20px;
        box-shadow: 0 10px 25px -5px rgba(0,0,0,0.15), 0 8px 10px -6px rgba(0,0,0,0.1);
        margin-bottom: 28px;
    }
    .hero-title {
        font-size: 2.25rem;
        font-weight: 800;
        letter-spacing: -0.03em;
        margin-bottom: 8px;
        background: linear-gradient(to right, #ffffff, #99f6e4);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    .hero-subtitle {
        color: #94a3b8;
        font-size: 1.05rem;
        font-weight: 400;
        line-height: 1.6;
    }

    /* Metric Cards */
    .metric-box {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 14px;
        padding: 20px 24px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        transition: transform 0.2s ease, box-shadow 0.2s ease;
    }
    .metric-box:hover {
        transform: translateY(-2px);
        box-shadow: 0 10px 15px -3px rgba(0,0,0,0.07);
    }
    .metric-label {
        color: #64748b;
        font-size: 0.85rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .metric-val {
        color: #0f172a;
        font-size: 1.85rem;
        font-weight: 700;
        margin-top: 6px;
    }
    .metric-sub {
        color: #10b981;
        font-size: 0.8rem;
        font-weight: 600;
        margin-top: 4px;
    }

    /* Diagnostic Status Banners */
    .diag-banner-normal {
        background: #ecfdf5;
        border-left: 6px solid #10b981;
        border-radius: 12px;
        padding: 22px 26px;
        margin: 20px 0;
    }
    .diag-banner-carrier {
        background: #fffbeb;
        border-left: 6px solid #f59e0b;
        border-radius: 12px;
        padding: 22px 26px;
        margin: 20px 0;
    }
    .diag-banner-pathogenic {
        background: #fef2f2;
        border-left: 6px solid #ef4444;
        border-radius: 12px;
        padding: 22px 26px;
        margin: 20px 0;
    }

    /* Professional Form Card */
    .auth-card {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 18px;
        padding: 34px 38px;
        box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05), 0 2px 4px -2px rgba(0,0,0,0.05);
    }

    /* Tabs Styling */
    .stTabs [data-baseweb="tab-list"] {
        border-bottom: 2px solid #e2e8f0;
        gap: 32px;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 12px 6px;
        font-size: 0.95rem;
        font-weight: 600;
        color: #64748b;
    }
    .stTabs [aria-selected="true"] {
        color: #0f766e !important;
        border-bottom-color: #0f766e !important;
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


def render_professional_auth():
    # Hero Title Bar
    st.markdown("""
    <div class="hero-banner">
        <div style="display:flex; justify-content:space-between; align-items:center;">
            <div>
                <div class="hero-title">GENOMIX · Clinical Motif Intelligence</div>
                <div class="hero-subtitle">
                    Linear-Time <b>Z-Algorithm O(N+M)</b> Sequence Telemetry & Machine Learning Diagnostic Portal.<br>
                    Integrated with Peer-Reviewed NCBI GenBank Triplet Repeat Expansion Registries.
                </div>
            </div>
            <div style="font-size: 3.5rem;">🧬</div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    auth_tab1, auth_tab2 = st.tabs(["🔐 Clinical Sign In", "📝 Create Practitioner Account"])

    with auth_tab1:
        col_c1, col_c2, col_c3 = st.columns([1, 2.2, 1])
        with col_c2:
            st.markdown('<div class="auth-card">', unsafe_allow_html=True)
            st.markdown("### Access Diagnostics Portal")
            st.caption("Enter your clinical or institutional researcher credentials to proceed.")

            with st.form("signin_form_pro"):
                login_id = st.text_input("Username or Institutional Email", placeholder="e.g. dr.smith@mayo.edu")
                login_pass = st.text_input("Password", type="password", placeholder="••••••••••••")
                remember = st.checkbox("Keep session active for 24 hours", value=True)
                submit_login = st.form_submit_button("Sign In to Diagnostic Console", type="primary", use_container_width=True)

                if submit_login:
                    if not login_id or not login_pass:
                        st.warning("Please specify both username and password.")
                    else:
                        user = authenticate_user(login_id, login_pass)
                        if user:
                            st.session_state.authenticated = True
                            st.session_state.user = user
                            st.toast(f"Welcome back, {user['full_name']}!", icon="✅")
                            st.rerun()
                        else:
                            st.error("Invalid credentials. Please verify your credentials or register below.")

            st.markdown('</div>', unsafe_allow_html=True)

    with auth_tab2:
        col_r1, col_r2, col_r3 = st.columns([1, 2.4, 1])
        with col_r2:
            st.markdown('<div class="auth-card">', unsafe_allow_html=True)
            st.markdown("### Register New Geneticist / Analyst")
            st.caption("Create an authenticated profile for genomic sequence tracking and automated auditing.")

            with st.form("signup_form_pro"):
                full_name = st.text_input("Full Name & Title", placeholder="e.g. Dr. Jane Smith, MD, FACMG")
                institution = st.text_input("Institution / Clinical Laboratory", placeholder="e.g. Department of Neurogenetics, Mayo Clinic")
                c_u1, c_u2 = st.columns(2)
                with c_u1:
                    new_user = st.text_input("Username", placeholder="e.g. jsmith266")
                with c_u2:
                    new_email = st.text_input("Work Email", placeholder="e.g. jsmith@hospital.org")
                new_pass = st.text_input("Security Passphrase", type="password", placeholder="Minimum 8 characters")
                agree = st.checkbox("I verify that all analysis adheres to patient genomic data protection guidelines.")

                submit_signup = st.form_submit_button("Register & Activate Credentials", use_container_width=True)

                if submit_signup:
                    if not agree:
                        st.warning("Please acknowledge data governance compliance.")
                    else:
                        ok, msg = register_user(new_user, new_email, new_pass, full_name, institution)
                        if ok:
                            st.success(f"{msg} You can now log in above.")
                        else:
                            st.error(msg)

            st.markdown('</div>', unsafe_allow_html=True)


def render_professional_dashboard():
    user = st.session_state.user

    # Sidebar Header
    with st.sidebar:
        st.markdown(f"""
        <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:12px; padding:16px; margin-bottom:18px;">
            <div style="font-size:0.75rem; font-weight:700; color:#64748b; text-transform:uppercase;">Active Clinician</div>
            <div style="font-size:1.15rem; font-weight:700; color:#0f172a; margin-top:3px;">{user['full_name']}</div>
            <div style="font-size:0.85rem; color:#0f766e; font-weight:500;">{user['institution'] if user['institution'] else 'Authorized Geneticist'}</div>
            <div style="font-size:0.75rem; color:#94a3b8; margin-top:4px;">ID: #{user['id']} · {user['email']}</div>
        </div>
        """, unsafe_allow_html=True)

        if st.button("🚪 Terminate Session (Sign Out)", use_container_width=True):
            st.session_state.authenticated = False
            st.session_state.user = None
            st.rerun()

        st.divider()
        st.markdown("### 🧬 Target Condition Profile")
        selected_disease_name = st.selectbox(
            "Select Clinical Panel",
            options=list(DISEASE_REGISTRY.keys()),
            format_func=lambda x: x.replace("_", " ")
        )
        disease_info = DISEASE_REGISTRY[selected_disease_name]

        st.markdown(f"""
        <div style="background:#ffffff; border:1px solid #cbd5e1; border-radius:10px; padding:14px; font-size:0.88rem;">
            <b>Gene Symbol:</b> <code>{disease_info['gene']}</code><br>
            <b>Target Motif:</b> <code>{disease_info['motif']}</code> (Triplet)<br>
            <b>Chromosome:</b> <code>{disease_info['chromosome']}</code><br>
            <b>Normal Threshold:</b> $\le {disease_info['clinical_thresholds']['normal_max']}$ repeats<br>
            <b>Pathogenic Limit:</b> $\ge {disease_info['clinical_thresholds']['pathogenic_min']}$ repeats<br>
            <div style="margin-top:8px; font-size:0.8rem; color:#64748b;">{disease_info['description']}</div>
        </div>
        """, unsafe_allow_html=True)

    # Top Banner
    st.markdown(f"""
    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #e2e8f0; padding-bottom:14px; margin-bottom:20px;">
        <div>
            <h2 style="margin:0; font-weight:800; color:#0f172a;">{selected_disease_name.replace('_', ' ')} Diagnostic Console</h2>
            <div style="color:#64748b; font-size:0.95rem;">Assay Engine: Linear-Time Z-Algorithm Verification & Multi-Feature Random Forest Classification</div>
        </div>
        <div style="background:#e0f2fe; color:#0369a1; padding:6px 14px; border-radius:20px; font-weight:600; font-size:0.85rem;">
            ACMG Standard Compliant
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Navigation Tabs
    tab_scan, tab_history, tab_benchmark, tab_methodology = st.tabs([
        "🔬 Sequence Diagnostic Scanner",
        "📈 User Analytics & Audit Trail",
        "🗄️ Peer-Reviewed NCBI Cohort",
        "📘 DAA Algorithmic Architecture"
    ])

    # TAB 1: DIAGNOSTIC SCANNER
    with tab_scan:
        source_mode = st.radio(
            "Select Sequence Input Modality:",
            [
                "📁 Pre-loaded NCBI Patient Clones",
                "🌐 Live NCBI GenBank Query (Entrez API)",
                "📤 Upload Custom FASTA File",
                "✍️ Direct Nucleotide Input"
            ],
            horizontal=True
        )

        sequence = ""
        sample_name = ""

        if source_mode == "📁 Pre-loaded NCBI Patient Clones":
            fasta_files = [f for f in os.listdir(DATA_DIR) if f.endswith(".fasta")]
            matched = [f for f in fasta_files if disease_info['gene'] in f or selected_disease_name in f]
            options = matched if matched else fasta_files
            selected_file = st.selectbox("Select Specimen Record", options)
            if selected_file:
                sample_name = selected_file
                header, sequence = read_fasta(os.path.join(DATA_DIR, selected_file))
                st.caption(f"**FASTA Record Header:** `{header}`")

        elif source_mode == "🌐 Live NCBI GenBank Query (Entrez API)":
            c_q1, c_q2 = st.columns([3, 1])
            with c_q1:
                acc_input = st.text_input("Enter NCBI Accession Number", value=disease_info["ncbi_accessions"][0])
            with c_q2:
                st.write("")
                st.write("")
                fetch_btn = st.button("Query NCBI GenBank", type="primary", use_container_width=True)

            if fetch_btn and acc_input:
                with st.spinner(f"Querying NCBI Entrez REST Service for '{acc_input}'..."):
                    try:
                        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc_input.strip()}&rettype=fasta&retmode=text"
                        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                        with urllib.request.urlopen(req) as resp:
                            fasta_raw = resp.read().decode('utf-8')
                        save_p = os.path.join(DATA_DIR, f"{acc_input.strip()}.fasta")
                        with open(save_p, "w", encoding="utf-8") as f:
                            f.write(fasta_raw)
                        header, sequence = read_fasta(save_p)
                        sample_name = f"NCBI Accession {acc_input.strip()}"
                        st.success(f"Verified & Loaded {len(sequence):,} base pairs from NCBI!")
                    except Exception as e:
                        st.error(f"NCBI Entrez Query Failure: {e}")

        elif source_mode == "📤 Upload Custom FASTA File":
            up_file = st.file_uploader("Select FASTA sequence file (.fasta, .fa)", type=["fasta", "fa"])
            if up_file:
                lines = up_file.read().decode("utf-8").splitlines()
                header = lines[0][1:] if lines[0].startswith(">") else "Local Sample"
                sequence = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
                sample_name = up_file.name

        else:
            raw_text = st.text_area("Paste DNA Nucleotides (A, C, G, T)", height=120, placeholder="ATGGCGACCCTGGAAAAGTCCTTCCAGCAGCAGCAGCAGCAGCAG...")
            if raw_text:
                sequence = "".join(raw_text.split()).upper()
                sample_name = "Direct Nucleotide Stream"

        # Action Button & Computation
        if sequence:
            st.divider()
            c_exec1, c_exec2 = st.columns([3, 1])
            with c_exec1:
                st.write(f"Specimen Ready: **{sample_name}** ({len(sequence):,} base pairs)")
            with c_exec2:
                run_calc = st.button("⚡ Execute Z-Scan & Predict", type="primary", use_container_width=True)

            if run_calc:
                with st.spinner("Processing Linear-Time Z-Algorithm [L, R] prefix sliding window..."):
                    motif = disease_info["motif"]
                    gene = disease_info["gene"]
                    thresh = disease_info["clinical_thresholds"]

                    # 1. Z-Algorithm Engine
                    feats = extract_z_features(sequence, motif, gene)
                    repeats = feats["max_consecutive_repeats"]

                    # 2. Machine Learning Classifier
                    clf = load_classifier()
                    df_feat = pd.DataFrame([feats])[FEATURE_COLUMNS]
                    pred_class = clf.predict(df_feat)[0]
                    pred_probs = clf.predict_proba(df_feat)[0]

                    labels = ["Normal (Healthy)", "Carrier (Intermediate)", "Reduced Penetrance", "Pathogenic (Full Penetrance)"]
                    predicted_label = labels[pred_class]
                    conf = pred_probs[pred_class] * 100.0

                    # 3. Clinical Classification
                    if repeats >= thresh["pathogenic_min"]:
                        clinical_status = "PATHOGENIC / FULL PENETRANCE"
                        risk_level = "High"
                        banner_class = "diag-banner-pathogenic"
                        icon = "🔴"
                    elif repeats > thresh["normal_max"]:
                        clinical_status = "INTERMEDIATE / CARRIER RISK"
                        risk_level = "Moderate"
                        banner_class = "diag-banner-carrier"
                        icon = "🟡"
                    else:
                        clinical_status = "NORMAL (WITHIN CLINICAL LIMITS)"
                        risk_level = "Low"
                        banner_class = "diag-banner-normal"
                        icon = "🟢"

                    # 4. Save to Database
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
                        "confidence_pct": round(conf, 2)
                    })

                    # Diagnostic Banner
                    st.markdown(f"""
                    <div class="{banner_class}">
                        <div style="font-size: 0.85rem; font-weight:700; color:#475569; text-transform:uppercase;">Official Clinical Assessment</div>
                        <div style="font-size: 1.6rem; font-weight:800; margin-top:4px;">{icon} {clinical_status} [{risk_level.upper()} RISK]</div>
                        <div style="font-size: 0.95rem; color:#334155; margin-top:6px;">
                            Identified <b>{repeats} consecutive '{motif}' repeats</b> spanning nucleotide loci <code>[{feats['locus_start']}:{feats['locus_end']}]</code>.
                            Machine Learning Classification: <b>{predicted_label}</b> with <b>{conf:.1f}%</b> confidence.
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

                    # Top KPI Metrics Cards
                    st.markdown("#### Primary Biomarkers & Algorithmic Telemetry")
                    m_c1, m_c2, m_c3, m_c4 = st.columns(4)

                    with m_c1:
                        st.markdown(f"""
                        <div class="metric-box">
                            <div class="metric-label">Max Consecutive Repeats</div>
                            <div class="metric-val">{repeats}</div>
                            <div class="metric-sub">Locus [{feats['locus_start']}:{feats['locus_end']}]</div>
                        </div>
                        """, unsafe_allow_html=True)

                    with m_c2:
                        st.markdown(f"""
                        <div class="metric-box">
                            <div class="metric-label">Total '{motif}' Matches</div>
                            <div class="metric-val">{feats['total_motif_occurrences']}</div>
                            <div class="metric-sub">{feats['motif_density_per_kb']} per kb</div>
                        </div>
                        """, unsafe_allow_html=True)

                    with m_c3:
                        st.markdown(f"""
                        <div class="metric-box">
                            <div class="metric-label">Checks Saved (Z-Box)</div>
                            <div class="metric-val">{feats['comparisons_saved']:,}</div>
                            <div class="metric-sub">{feats['efficiency_ratio']*100:.1f}% vs Naive Scan</div>
                        </div>
                        """, unsafe_allow_html=True)

                    with m_c4:
                        st.markdown(f"""
                        <div class="metric-box">
                            <div class="metric-label">Z-Score Enrichment</div>
                            <div class="metric-val">{feats['z_enrichment_score']} σ</div>
                            <div class="metric-sub">Above Random Background</div>
                        </div>
                        """, unsafe_allow_html=True)

                    st.write("")

                    # Comparative Visuals (Plotly)
                    c_viz1, c_viz2 = st.columns(2)

                    with c_viz1:
                        st.markdown("##### ⚡ DAA Prefix Window Optimization")
                        fig_opt = go.Figure(data=[
                            go.Bar(name='Direct Checks Executed', x=['Z-Algorithm'], y=[feats['z_comparisons']], marker_color='#0f766e'),
                            go.Bar(name='Redundant Checks Saved', x=['Z-Algorithm'], y=[feats['comparisons_saved']], marker_color='#99f6e4')
                        ])
                        fig_opt.update_layout(barmode='stack', height=300, margin=dict(l=20, r=20, t=30, b=20),
                                              legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
                        st.plotly_chart(fig_opt, use_container_width=True)

                    with c_viz2:
                        st.markdown("##### 🧠 Multi-Class Model Probability Distribution")
                        df_probs = pd.DataFrame({"Clinical State": labels, "Probability (%)": pred_probs * 100.0})
                        fig_p = px.bar(
                            df_probs, x="Probability (%)", y="Clinical State", orientation='h',
                            color="Probability (%)", color_continuous_scale="Teal", text="Probability (%)", range_x=[0, 100]
                        )
                        fig_p.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
                        fig_p.update_layout(height=300, margin=dict(l=20, r=20, t=30, b=20), coloraxis_showscale=False)
                        st.plotly_chart(fig_p, use_container_width=True)

    # TAB 2: USER AUDIT & ANALYTICS
    with tab_history:
        st.header("Clinician Activity Dashboard & Patient Analytics")
        history = get_user_analyses(user["id"])

        if not history:
            st.info("No recorded diagnostics in your profile yet. Run your first sequence in the scanner to populate telemetry.")
        else:
            df_h = pd.DataFrame(history)

            # High-level Analytics KPIs
            k1, k2, k3, k4 = st.columns(4)
            k1.metric("Screenings Completed", len(df_h))
            k2.metric("Positive Detections", len(df_h[df_h["risk_level"] == "High"]))
            k3.metric("Total Character Skips", f"{df_h['comparisons_saved'].sum():,}")
            k4.metric("Avg Repeat Length", f"{df_h['max_repeats'].mean():.1f} repeats")

            st.divider()

            c_chart1, c_chart2 = st.columns(2)

            with c_chart1:
                st.subheader("Patient Cohort Risk Stratification")
                risk_df = df_h["risk_level"].value_counts().reset_index()
                risk_df.columns = ["Risk", "Count"]
                fig_pie = px.donut = px.pie(
                    risk_df, names="Risk", values="Count", hole=0.5,
                    color="Risk", color_discrete_map={"Low": "#10b981", "Moderate": "#f59e0b", "High": "#ef4444"}
                )
                fig_pie.update_layout(height=320, margin=dict(l=10, r=10, t=10, b=10))
                st.plotly_chart(fig_pie, use_container_width=True)

            with c_chart2:
                st.subheader("Repeat Expansions Across Verified Specimens")
                fig_hist_bar = px.bar(
                    df_h, x="sample_name", y="max_repeats", color="risk_level",
                    color_discrete_map={"Low": "#10b981", "Moderate": "#f59e0b", "High": "#ef4444"},
                    labels={"sample_name": "Specimen ID", "max_repeats": "Consecutive Repeats"}
                )
                fig_hist_bar.update_layout(height=320, margin=dict(l=10, r=10, t=10, b=10))
                st.plotly_chart(fig_hist_bar, use_container_width=True)

            st.subheader("Comprehensive Diagnostic Audit Trail")
            st.dataframe(
                df_h[[
                    "id", "sample_name", "disease", "target_motif",
                    "seq_length", "max_repeats", "prediction_label",
                    "risk_level", "confidence_pct", "created_at"
                ]],
                use_container_width=True
            )

    # TAB 3: REAL NCBI COHORT BENCHMARK
    with tab_benchmark:
        st.header("Peer-Reviewed NCBI GenBank Triplet Repeat Cohort")
        st.caption("All records below are 100% genuine biological specimens fetched from NCBI GenBank.")

        summary_file = os.path.join(RESULTS_DIR, "real_ncbi_analysis_summary.json")
        if os.path.exists(summary_file):
            with open(summary_file, "r", encoding="utf-8") as jf:
                s_data = json.load(jf)

            table_rows = []
            for item in s_data:
                f_item = item["features"]
                table_rows.append({
                    "NCBI Accession": item["accession"],
                    "Target Gene": item["gene"],
                    "Condition": item["disease"].replace("_", " "),
                    "Motif": item["target_motif"],
                    "Length (bp)": f"{item['sequence_length_bp']:,}",
                    "Max Tandem Repeats": f_item["max_consecutive_repeats"],
                    "Z-Comparisons": f"{f_item['z_comparisons']:,}",
                    "Checks Bypassed": f"{f_item['comparisons_saved']:,}",
                    "Efficiency": f"{f_item['efficiency_ratio']*100:.1f}%",
                    "Diagnosis": item["clinical_status"]
                })

            st.dataframe(pd.DataFrame(table_rows), use_container_width=True)

    # TAB 4: METHODOLOGY
    with tab_methodology:
        st.header("Algorithmic & Theoretical Foundations")
        st.markdown(r"""
        ### 1. The Classical Z-Algorithm
        For any input string $S$ of length $N$, the algorithm computes the fundamental array $Z$ in **strictly linear time $O(N)$**:
        $$Z[i] = \max \{ k \mid S[i \dots i+k-1] = S[0 \dots k-1] \}$$

        To locate exact occurrences of motif $P$ in genomic sequence $T$, we formulate the concatenated string:
        $$S = P \mathbin{\$} T$$
        where $\$$ is a unique delimiter not present in the alphabet $\{A, C, G, T\}$.

        ### 2. Algorithmic Prefix Window $[L, R]$
        Instead of re-comparing characters from scratch at each position, the Z-algorithm maintains the rightmost match segment $[L, R]$:
        - If the current index $i \le R$, we compute $k = i - L$.
        - If $Z[k] < R - i + 1$, the value $Z[i]$ is copied in **$O(1)$ time with zero character comparisons**.
        - Only when $Z[k] \ge R - i + 1$ does the algorithm perform explicit character checks, extending the right bound $R$.
        """)


def main():
    if not st.session_state.authenticated:
        render_professional_auth()
    else:
        render_professional_dashboard()


if __name__ == "__main__":
    main()
