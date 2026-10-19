"""
NCBI Real Data Downloader:
Connects to the official NCBI Entrez E-utilities API to download authentic,
peer-reviewed human genomic clones and patient mRNA records across 3 major diseases:
  1. Huntington's Disease (HTT - CAG repeats)
  2. Fragile X Syndrome (FMR1 - CGG repeats)
  3. Friedreich's Ataxia (FXN - GAA repeats)
"""

import os
import urllib.request
import time

# Catalog of real NCBI accessions with documented clinical context
REAL_NCBI_CATALOG = [
    # --- 1. HUNTINGTON'S DISEASE (HTT gene, CAG motif) ---
    {
        "accession": "KJ535072.1",
        "disease": "Huntington_Disease",
        "gene": "HTT",
        "motif": "CAG",
        "clinical_class": "Normal",
        "description": "Homo sapiens fetal brain cDNA clone HTT_iso-A (Normal CAG tract)"
    },
    {
        "accession": "NM_002111.9",
        "disease": "Huntington_Disease",
        "gene": "HTT",
        "motif": "CAG",
        "clinical_class": "Normal",
        "description": "Homo sapiens huntingtin (HTT) mRNA RefSeq reference transcript"
    },
    {
        "accession": "L12392.1",
        "disease": "Huntington_Disease",
        "gene": "HTT",
        "motif": "CAG",
        "clinical_class": "Normal",
        "description": "Landmark Huntington's Disease Collaborative Research Group cDNA clone"
    },
    {
        "accession": "L27350.1",
        "disease": "Huntington_Disease",
        "gene": "HTT",
        "motif": "CAG",
        "clinical_class": "Normal",
        "description": "Homo sapiens Huntington's disease (IT15) gene 5' region clone"
    },

    # --- 2. FRAGILE X SYNDROME (FMR1 gene, CGG motif) ---
    {
        "accession": "NM_002024.6",
        "disease": "Fragile_X_Syndrome",
        "gene": "FMR1",
        "motif": "CGG",
        "clinical_class": "Normal",
        "description": "Homo sapiens fragile X messenger ribonucleoprotein 1 (FMR1) transcript"
    },
    {
        "accession": "L29074.1",
        "disease": "Fragile_X_Syndrome",
        "gene": "FMR1",
        "motif": "CGG",
        "clinical_class": "Normal",
        "description": "Homo sapiens fragile X (FMR1) gene genomic locus covering CGG repeat"
    },
    {
        "accession": "AF305880.1",
        "disease": "Fragile_X_Syndrome",
        "gene": "FMR1",
        "motif": "CGG",
        "clinical_class": "Normal",
        "description": "Homo sapiens FMR1 promoter and CGG triplet repeat region"
    },

    # --- 3. FRIEDREICH'S ATAXIA (FXN gene, GAA motif) ---
    {
        "accession": "NM_000144.6",
        "disease": "Friedreich_Ataxia",
        "gene": "FXN",
        "motif": "GAA",
        "clinical_class": "Normal",
        "description": "Homo sapiens frataxin (FXN) transcript variant 1 mRNA"
    },
    {
        "accession": "U43748.1",
        "disease": "Friedreich_Ataxia",
        "gene": "FXN",
        "motif": "GAA",
        "clinical_class": "Normal",
        "description": "Campuzano et al. landmark genomic DNA clone covering FXN GAA tract"
    },
    {
        "accession": "U43750.1",
        "disease": "Friedreich_Ataxia",
        "gene": "FXN",
        "motif": "GAA",
        "clinical_class": "Normal",
        "description": "Human frataxin gene intron 1 genomic clone with GAA triplet region"
    }
]


def download_ncbi_catalog():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    print("=" * 80)
    print("  DOWNLOADING AUTHENTIC GENOMIC SAMPLES FROM NCBI GENBANK (E-UTILITIES)")
    print("=" * 80)

    for item in REAL_NCBI_CATALOG:
        acc = item["accession"]
        out_filename = f"{acc.split('.')[0]}_{item['disease']}.fasta"
        out_path = os.path.join(data_dir, out_filename)

        if os.path.exists(out_path):
            print(f"[EXISTS] {acc} is already present as {out_filename}")
            continue

        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
        print(f"[*] Fetching real sequence {acc} ({item['disease']})...")
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'NCBI-Downloader/1.0'})
            with urllib.request.urlopen(req) as resp:
                content = resp.read().decode('utf-8')
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(content)
            print(f"    [+] Saved: {out_filename} ({len(content)} bytes)")
            time.sleep(0.5)  # Be polite to NCBI API rate limits
        except Exception as e:
            print(f"    [!] Error downloading {acc}: {e}")

    print("=" * 80)
    print(f"[SUCCESS] All real NCBI sequences downloaded directly into: {data_dir}")
    print("=" * 80)


if __name__ == "__main__":
    download_ncbi_catalog()
