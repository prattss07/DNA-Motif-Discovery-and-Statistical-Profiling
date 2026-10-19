"""
Multi-Disease Genetic Catalog:
Contains verified NCBI GenBank accessions and repeat motifs for major
genetic diseases, allowing dynamic addition of new samples at evaluation time.
"""

DISEASE_REGISTRY = {
    "Huntington_Disease": {
        "gene": "HTT",
        "motif": "CAG",
        "chromosome": "4p16.3",
        "clinical_thresholds": {
            "normal_max": 26,
            "intermediate_max": 35,
            "reduced_penetrance_max": 39,
            "pathogenic_min": 40
        },
        "description": "Progressive neurodegenerative chorea caused by polyglutamine expansion in HTT gene.",
        "ncbi_accessions": [
            "NM_002111.9",   # Human reference HTT mRNA
            "KJ535072.1",    # Brain cDNA clone
            "L12392.1",      # Landmark HD study clone
            "AH002939.2"     # Exon 1 polymorphic region
        ]
    },
    "Fragile_X_Syndrome": {
        "gene": "FMR1",
        "motif": "CGG",
        "chromosome": "Xq27.3",
        "clinical_thresholds": {
            "normal_max": 44,
            "intermediate_max": 54,
            "reduced_penetrance_max": 199,
            "pathogenic_min": 200
        },
        "description": "Leading inherited cause of intellectual disability and autism spectrum traits.",
        "ncbi_accessions": [
            "NM_002024.6",   # Human FMR1 transcript variant 1
            "L29074.1",      # FMR1 5' untranslated CGG repeat region
            "AF305880.1"     # Exon 1 CGG repeat region
        ]
    },
    "Friedreich_Ataxia": {
        "gene": "FXN",
        "motif": "GAA",
        "chromosome": "9q21.11",
        "clinical_thresholds": {
            "normal_max": 33,
            "intermediate_max": 65,
            "reduced_penetrance_max": 99,
            "pathogenic_min": 100
        },
        "description": "Autosomal recessive spinocerebellar ataxia causing sensory loss and cardiomyopathy.",
        "ncbi_accessions": [
            "NM_000144.5",   # Human FXN frataxin transcript
            "U43750.1",      # Frataxin intron 1 GAA triplet repeat region
            "AF048560.1"     # GAA repeat polymorphic allele
        ]
    },
    "Myotonic_Dystrophy_Type_1": {
        "gene": "DMPK",
        "motif": "CTG",
        "chromosome": "19q13.32",
        "clinical_thresholds": {
            "normal_max": 34,
            "intermediate_max": 49,
            "reduced_penetrance_max": 99,
            "pathogenic_min": 100
        },
        "description": "Multisystem disorder characterized by myotonia, progressive muscle weakness and cataracts.",
        "ncbi_accessions": [
            "NM_001081560.3", # Human DMPK transcript
            "L08835.1",       # DMPK 3' untranslated CTG region
            "M87338.1"        # Myotonin protein kinase CTG repeats
        ]
    },
    "Spinocerebellar_Ataxia_Type_1": {
        "gene": "ATXN1",
        "motif": "CAG",
        "chromosome": "6p22.3",
        "clinical_thresholds": {
            "normal_max": 35,
            "intermediate_max": 38,
            "reduced_penetrance_max": 43,
            "pathogenic_min": 44
        },
        "description": "Progressive neurodegenerative ataxia affecting balance, speech and coordination.",
        "ncbi_accessions": [
            "NM_000332.4",   # Human ATXN1 ataxin 1 transcript
            "L16628.1"       # Ataxin-1 CAG repeat region
        ]
    }
}
