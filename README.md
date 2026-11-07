# DNA Motif Discovery & Multi-Disease Genetic Profiling Platform (C++ Edition)

A high-performance computational biology and clinical diagnostic suite implemented in **pure modern C++ (C++17)** utilizing the fundamental **linear-time Z-Algorithm ($O(N + M)$)** for exact DNA motif identification, tandem repeat expansion profiling, statistical background modeling, and clinical risk prediction across major human genetic disorders.

Built with 100% authentic biological sequences sourced from **NCBI GenBank**, this platform includes an embedded high-performance **HTTP Web Server & Interactive Dashboard** with practitioner authentication (Sign Up, Sign In, Logout), real-time repeat telemetry, and persistent user analysis audit logs—**entirely in C++ with zero external runtime dependencies**.

---

## ⚡ Key Highlights of the C++ Implementation

* **Zero External Dependencies:** Native C++17 implementation with custom WinSock2 HTTP server, standard threading, and zero Python/pip dependencies.
* **Deterministic $O(N + M)$ Z-Algorithm:** Implements the classical sliding prefix window $[L, R]$ to eliminate redundant character comparisons, achieving over 55% comparison efficiency over naive scanning.
* **ACMG Diagnostic Registry:** Automated clinical categorization into **Normal (Healthy)**, **Carrier (Intermediate)**, **Reduced Penetrance**, and **Pathogenic (Full Penetrance)** across 5 major human triplet repeat expansion disorders.
* **Dual Operation Modes:**
  1. **Web Dashboard Mode:** Embedded browser-based interface (`http://localhost:8080`) with authentication, file selection, direct DNA pasting, and structured pathology reporting.
  2. **High-Speed CLI Mode:** Single-file screening, batch benchmarking across all 10 NCBI GenBank reference clones, and interactive console menus.

---

## 🧬 Targeted Genetic Conditions & Biomarkers

The system profiles polymorphic trinucleotide repeat expansions linked to neurological and multisystem diseases according to American College of Medical Genetics (ACMG) diagnostic standards:

| Disease | Gene | Chromosome | Target Motif | Normal Range | Pathogenic Threshold |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Huntington's Disease (HD)** | *HTT* | 4p16.3 | `CAG` | $\le 26$ repeats | $\ge 40$ repeats |
| **Fragile X Syndrome (FXS)** | *FMR1* | Xq27.3 | `CGG` | $\le 44$ repeats | $\ge 200$ repeats |
| **Friedreich's Ataxia (FRDA)** | *FXN* | 9q21.11 | `GAA` | $\le 33$ repeats | $\ge 100$ repeats |
| **Myotonic Dystrophy Type 1** | *DMPK* | 19q13.32 | `CTG` | $\le 34$ repeats | $\ge 100$ repeats |
| **Spinocerebellar Ataxia 1** | *ATXN1* | 6p22.3 | `CAG` | $\le 35$ repeats | $\ge 44$ repeats |

---

## 📁 Authentic NCBI GenBank Dataset Cohort

All sequences analyzed are genuine, peer-reviewed records stored in `data/`:

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

## 🔬 Core Algorithmic Architecture

### 1. The Classical Linear-Time Z-Algorithm ($O(N + M)$)
For a concatenated string $S = P \mathbin{\$} T$ of length $|P| + 1 + |T|$ (where $P$ is the motif and $T$ is the genome sequence), the algorithm computes the fundamental array $Z$:
$$Z[i] = \max \{ k \mid S[i \dots i+k-1] = S[0 \dots k-1] \}$$

* **Active Window $[L, R]$:** Maintains the rightmost match segment. When current index $i \le R$, previous computations are reused in $O(1)$ time, skipping redundant character comparisons.
* **Tandem Repeat Clustering:** Discovered match indices are evaluated for consecutive adjacency ($k, k+3, k+6, \dots$) to detect unbroken repeat chains.


### 2. Statistical Background Model & Risk Prediction
* **Markov Null Model:** Calculates expected motif counts based on nucleotide background composition ($GC$-ratio).
* **Z-Score Enrichment:** Measures statistical deviation of observed motif hits above random expectation:
  $$Z_{\text{enrichment}} = \frac{\text{Observed} - E}{\sigma}$$
* **Probabilistic Classification:** Multi-feature classifier calibrated across all four risk classes (**Normal**, **Carrier**, **Reduced Penetrance**, **Pathogenic**) with confidence percentages and actionable clinical recommendations.

---

## 📁 Project Structure

```
D:\PACKAGES\SEM5\DAA\
│
├── include/                            # C++ Header declarations
│   ├── z_algorithm.hpp                 # Linear-time Z-Algorithm O(N+M)
│   ├── fasta_parser.hpp                # Robust multi-line FASTA parser
│   ├── tandem_analyzer.hpp             # Adjacent repeat expansion clustering
│   ├── statistical_model.hpp           # Markov Null model & Z-enrichment
│   ├── disease_registry.hpp            # ACMG clinical diagnostic criteria
│   ├── classifier.hpp                  # Probabilistic clinical classification
│   ├── auth_system.hpp                 # Pure C++ SHA-256 & user session store
│   ├── pathology_report.hpp            # Structured clinical report formatter
│   └── web_server.hpp                  # Embedded WinSock2 HTTP server
│
├── cpp_src/                            # C++ Source implementations
│   ├── z_algorithm.cpp
│   ├── fasta_parser.cpp
│   ├── tandem_analyzer.cpp
│   ├── statistical_model.cpp
│   ├── disease_registry.cpp
│   ├── classifier.cpp
│   ├── auth_system.cpp
│   ├── pathology_report.cpp
│   ├── web_server.cpp
│   └── main.cpp                        # Unified CLI & server entry point
│
├── web/                                # Web UI frontend
│   └── index.html                      # Modern, clean, dark-mode clinical dashboard
│
├── data/                               # 100% Authentic NCBI GenBank FASTA sequences
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
├── bin/
│   └── genomix.exe                     # Statically linked standalone binary
│
├── build.bat                           # 1-Click build script (GCC MinGW-W64)
├── Makefile                            # Standard Make build file
└── README.md                           # Documentation
```

---


.\bin\genomix.exe
```

---


