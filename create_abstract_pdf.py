"""
Generate DAA Package Abstract PDF matching the exact font, size, alignment and layout
of the original scanned document.
"""

import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_RIGHT, TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib import colors

def generate_pdf():
    pdf_path = os.path.join(r"D:\PACKAGES\SEM5\DAA", "DAA_PACKAGE_ABSTRACT.pdf")
    
    # 1 inch = 72 points
    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=A4,
        leftMargin=65,
        rightMargin=65,
        topMargin=70,
        bottomMargin=70
    )
    
    story = []
    
    # Header Style
    header_style = ParagraphStyle(
        name='AbstractHeader',
        fontName='Times-Bold',
        fontSize=16,
        leading=22,
        alignment=TA_CENTER,
        textColor=colors.black,
        spaceAfter=26
    )
    
    # Package Title Style
    title_style = ParagraphStyle(
        name='PackageTitle',
        fontName='Times-Bold',
        fontSize=11,
        leading=16,
        alignment=TA_LEFT,
        textColor=colors.black,
        spaceAfter=20
    )
    
    # Body Paragraph Style
    body_style = ParagraphStyle(
        name='AbstractBody',
        fontName='Times-Roman',
        fontSize=11,
        leading=15.5,
        alignment=TA_JUSTIFY,
        textColor=colors.black,
        spaceAfter=14
    )
    
    # Team Members Header Style
    team_title_style = ParagraphStyle(
        name='TeamTitle',
        fontName='Times-Bold',
        fontSize=11.5,
        leading=16,
        alignment=TA_RIGHT,
        textColor=colors.black,
        spaceAfter=3
    )
    
    # Team Member Name Style
    team_name_style = ParagraphStyle(
        name='TeamMember',
        fontName='Times-Roman',
        fontSize=10.5,
        leading=15,
        alignment=TA_RIGHT,
        textColor=colors.black
    )
    
    # 1. Header
    story.append(Paragraph("<b>DAA PACKAGE ABSTRACT</b>", header_style))
    
    # 2. Package Title
    story.append(Paragraph("<b><i>PACKAGE TITLE:</i> DNA MOTIF DISCOVERY AND STATISTICAL PROFILING</b>", title_style))
    
    # 3. Paragraph 1
    p1 = (
        "The DNA Motif Discovery and Statistical Profiling system is a Design and Analysis of Algorithms "
        "(DAA) project designed to demonstrate the efficient identification and analysis of exact nucleotide "
        "patterns in genomic sequences using the linear-time Z-Algorithm. The project focuses on discovering "
        "biologically significant DNA motifs across real genome samples and evaluating their frequency distributions "
        "to establish their role in molecular biology and disease diagnosis."
    )
    story.append(Paragraph(p1, body_style))
    
    # 4. Paragraph 2
    p2 = (
        "The application models the linear-time string matching process by constructing a concatenated string of the "
        "motif pattern, a unique sentinel delimiter, and the genomic text. It illustrates the fundamental working of the "
        "Z-Algorithm, which computes the length of the longest substring starting at each position that matches a prefix "
        "of the string in <i>O(n + m)</i> time. By maintaining an active matching window defined by [L, R] boundaries, the "
        "algorithm reuses previously computed prefix matches to avoid redundant character comparisons, dramatically "
        "outperforming naive search methods on highly repetitive genomic sequences."
    )
    story.append(Paragraph(p2, body_style))
    
    # 5. Paragraph 3
    p3 = (
        "The project aims to provide an interactive and visual representation of motif occurrences and statistical profiles "
        "across diverse biological sequences sourced from NCBI GenBank. It incorporates a Markov null model calibrated to "
        "background nucleotide (GC) content to compute statistical enrichment Z-scores and identify significant tandem repeat "
        "expansions. Users can observe algorithmic telemetry, motif distributions, and comparative profiles across different "
        "specimens to classify samples into normal, carrier, and pathogenic risk categories."
    )
    story.append(Paragraph(p3, body_style))
    
    # 6. Paragraph 4
    p4 = (
        "This project helps develop a practical understanding of advanced linear-time string matching algorithms, algorithm "
        "efficiency optimization, and statistical data modeling in computational biology. It can serve as an educational and "
        "analytical tool for understanding how fundamental algorithmic techniques like the Z-Algorithm can be leveraged to "
        "analyze massive genomic datasets and predict genetic disorders reliably."
    )
    story.append(Paragraph(p4, body_style))
    
    # 7. Team Members Section
    story.append(Spacer(1, 35))
    story.append(Paragraph("<b>TEAM MEMBERS</b>", team_title_style))
    story.append(Paragraph("PRATHAMESH P(24PD28)", team_name_style))
    
    doc.build(story)
    print(f"SUCCESS: Generated PDF at {pdf_path}")

if __name__ == "__main__":
    generate_pdf()
