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





