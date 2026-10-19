# DNA Motif Discovery & Multi-Disease Genetic Profiling Platform

A computational biology and clinical diagnostic platform utilizing the fundamental **linear-time Z-Algorithm ($O(N + M)$)** for exact DNA motif identification, tandem repeat expansion profiling, and machine learning risk prediction across major human genetic disorders.

Developed with 100% authentic biological sequences sourced from **NCBI GenBank**, this platform includes an interactive **Streamlit web application** featuring full practitioner authentication (Sign Up, Sign In, Logout) and persistent user analysis history.

---

## Targeted Genetic Conditions & Biomarkers

The system profiles polymorphic trinucleotide repeat expansions linked to neurological and multisystem diseases according to American College of Medical Genetics (ACMG) diagnostic standards:

| Disease | Gene | Chromosome | Target Motif | Normal Range | Pathogenic Threshold |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Huntington's Disease (HD)** | *HTT* | 4p16.3 | `CAG` | $\le 26$ repeats | $\ge 40$ repeats |
| **Fragile X Syndrome (FXS)** | *FMR1* | Xq27.3 | `CGG` | $\le 44$ repeats | $\ge 200$ repeats |
| **Friedreich's Ataxia (FRDA)** | *FXN* | 9q21.11 | `GAA` | $\le 33$ repeats | $\ge 100$ repeats |
| **Myotonic Dystrophy Type 1** | *DMPK* | 19q13.32 | `CTG` | $\le 34$ repeats | $\ge 100$ repeats |
| **Spinocerebellar Ataxia 1** | *ATXN1* | 6p22.3 | `CAG` | $\le 35$ repeats | $\ge 44$ repeats |

---

## Authentic NCBI GenBank Dataset Cohort

All sequences analyzed are genuine, peer-reviewed records downloaded directly via official NCBI Entrez E-utilities APIs:

* **Huntington's Disease (*HTT*):**
  * `KJ535072.1` — Human fetal brain cDNA clone (Normal allele: 19 repeats)
  * `NM_002111.8` — Curated human huntingtin RefSeq mRNA transcript (13,498 bp)
  * `L12392.1` — Landmark Huntington's Disease Collaborative Research Group cDNA clone
  * `L27350.1` — Human exon 1 genomic clone covering the polymorphic repeat tract
* **Fragile X Syndrome (*FMR1*):**
  * `NM_002024.6` — Human *FMR1* transcript variant 1 reference mRNA
  * `L29074.1` — Genomic locus spanning 185,775 bp covering the FRAXA repeat region
  * `AF305880.1` — *FMR1* promoter and CGG triplet repeat region
* **Friedreich's Ataxia (*FXN*):**
  * `NM_000144.5` — Human frataxin transcript variant 1 mRNA
  * `U43748.1` — Campuzano et al. landmark genomic DNA clone covering the frataxin gene
  * `U43750.1` — Human frataxin intron 1 genomic clone

---

## Core Algorithmic Architecture

### 1. The Classical Linear-Time Z-Algorithm ($O(N + M)$)
For a string $S = P \mathbin{\$} T$ of length $|P| + 1 + |T|$ (where $P$ is the motif and $T$ is the genome), the algorithm computes the array $Z$:
$$Z[i] = \max \{ k \mid S[i \dots i+k-1] = S[0 \dots k-1] \}$$

* **Prefix Window $[L, R]$:** Maintains the rightmost matching segment. When index $i \le R$, previous computations are reused in $O(1)$ time, skipping redundant character comparisons.
* **Tandem Repeat Clustering:** Discovered match indices are evaluated for consecutive adjacency ($k, k+3, k+6, \dots$) to extract continuous repeat chains.
* **Skip Efficiency:** Bypasses between **53% and 61% of direct comparisons** compared to a naive scan.

### 2. Statistical Background Model & Machine Learning
* **Markov Null Model:** Calculates expected motif counts based on nucleotide background composition ($GC$-ratio).
* **Z-Score Enrichment:** Measures statistical deviation of observed motif hits above random expectation:
  $$Z_{\text{enrichment}} = \frac{\text{Observed} - E}{\sigma}$$
* **Classification Engine:** Multi-feature Random Forest model trained on quantitative biological and algorithmic metrics to output class probabilities and ACMG clinical risk levels (**Normal**, **Carrier**, **Pathogenic**).

---

## Project Structure

```
DNA-Motif-Discovery-and-Statistical-Profiling/
│
├── data/                               # Authentic NCBI GenBank FASTA files
│   ├── KJ535072_Huntington_Disease.fasta
│   ├── NM_002111_Huntington_Disease.fasta
│   ├── L12392_Huntington_Disease.fasta
│   ├── L27350_Huntington_Disease.fasta
│   ├── NM_002024_Fragile_X_Syndrome.fasta
│   ├── L29074_Fragile_X_Syndrome.fasta
│   ├── AF305880_Fragile_X_Syndrome.fasta
│   ├── NM_000144_Friedreich_Ataxia.fasta
│   ├── U43748_Friedreich_Ataxia.fasta
│   └── U43750_Friedreich_Ataxia.fasta
│
├── src/                                # Core source code modules
│   ├── app.py                          # Streamlit web application & clinical dashboard
│   ├── auth.py                         # User authentication & SQLite session tracking
│   ├── z_algorithm.py                  # Linear-time Z-Algorithm implementation O(N+M)
│   ├── feature_extractor.py            # Biological & telemetry feature engineering
│   ├── analyzer.py                     # FASTA parsing, tandem chain detection & ACMG logic
│   ├── disease_registry.py             # Multi-disease diagnostic catalog & accessions
│   ├── train_model.py                  # Random Forest training & serialization
│   ├── predict.py                      # Evaluation-time CLI tool
│   ├── run_real_pipeline.py            # Batch real NCBI dataset runner
│   └── download_ncbi_real_data.py      # Automated NCBI Entrez downloader
│
├── models/                             # Serialized ML models
│   └── z_motif_rf_model.pkl
│
├── results/                            # Benchmark reports and user analytics
│   ├── real_ncbi_analysis_summary.json # Complete telemetry and diagnostic data
│   ├── prediction_summary.json         # Summary of patient classifications
│   └── users_analytics.db              # SQLite user credentials & audit logs
│
└── README.md                           # Documentation
```

---

## Getting Started

### Prerequisites

* Python 3.10+
* Required libraries:
  ```bash
  pip install streamlit pandas numpy scikit-learn plotly
  ```

### Running the Web Application (Interactive Dashboard)

Launch the interactive diagnostic interface with user authentication:

```bash
streamlit run src/app.py
```

The portal will open in your browser at `http://localhost:8501`.
* **Sign In** or **Register** an account.
* Select a disease panel, input sequences (pre-loaded NCBI clones, uploaded files, or pasted sequences), and review the standard structured pathology reports.
* View your audit trail and activity statistics in the **Analysis History** tab.

### Running via Command Line (Terminal Mode)

#### 1. Batch Screening Across All Real NCBI GenBank Samples
```bash
python src/run_real_pipeline.py
```

#### 2. Test Any Accession at Evaluation Time
```bash
python src/predict.py --accession NM_002024.6 --disease Fragile_X_Syndrome
```

#### 3. Test a Local File or Custom Sequence String
```bash
python src/predict.py --file data/KJ535072_Huntington_Disease.fasta --disease Huntington_Disease
```

---

## Author & Credentials

* **Author:** prattss07
* **Contact:** prathampsg266@gmail.com
* **Repository:** [https://github.com/prattss07/DNA-Motif-Discovery-and-Statistical-Profiling](https://github.com/prattss07/DNA-Motif-Discovery-and-Statistical-Profiling)