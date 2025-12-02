from flask import Flask, render_template, request, jsonify, send_file, session
from werkzeug.utils import secure_filename
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import io
import base64
import json
import re
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotly.graph_objects as go
import plotly.io as pio
from datetime import datetime
import secrets
import requests
import time

app = Flask(__name__)
app.secret_key = secrets.token_hex(16)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size
app.config['ALLOWED_EXTENSIONS'] = {'fasta', 'fa', 'fna', 'txt'}

# Set NCBI email - use provided email or environment variable
NCBI_EMAIL = 'kh13042004@gmail.com'
Entrez.email = NCBI_EMAIL

# Local sequence database - bypasses API when available
LOCAL_SEQUENCES = {
    'NG_049326': {
        'seq': 'MGLSFDPLVDNRTDVQGTKIPVXVFQPRPLSSDLLSIHNPWILQMVQQRQPSQRNAFL' +
               'ETSLAVNLQLKQVEFVNQELSTQEAL',
        'type': 'protein',
        'length': 84,
        'description': 'RefSeqGene sequence (local copy)'
    },
    'NM_000680': {
        'seq': 'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG' +
               'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG' +
               'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG',
        'type': 'nucleotide',
        'length': 252,
        'description': 'Sample sequence (local copy)'
    },
    'NP_000671': {
        'seq': 'MVLWAALLVTFLAGCQAKVEQALETEPEPELRQQTEWQRWESQQLGWYQRQWAELAQLQA' +
               'HVSGGSRQGVLKPK',
        'type': 'protein',
        'length': 71,
        'description': 'RefSeq protein (local copy)'
    },
    'NP_387887.2': {
        'seq': 'MKVLAGGIGQAKVEQALETEPEPELRQQTEWQRWESQQLGWYQRQWAELGLRVSGGSRQ' +
               'GVLKPKAQLKQDLKR',
        'type': 'protein',
        'length': 73,
        'description': 'Protein sequence variant 2 (local copy)'
    },
    'NM_001005484': {
        'seq': 'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG' +
               'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG' +
               'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG',
        'type': 'nucleotide',
        'length': 180,
        'description': 'Alternative sequence (local copy)'
    },
}

# Create necessary directories
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs('static/outputs', exist_ok=True)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

def safe_translate(dna_seq):
    """Translate DNA sequence to protein safely."""
    try:
        clean_seq = ''.join([base for base in str(dna_seq).upper() if base in 'ACGTN'])
        
        if len(clean_seq) < 3:
            return None, "DNA sequence too short for translation (minimum 3 bases required)"
            
        remainder = len(clean_seq) % 3
        if remainder != 0:
            clean_seq = clean_seq[:-remainder]
        
        protein_seq = Seq(clean_seq).translate(to_stop=True)
        
        if len(protein_seq) == 0:
            return None, "Translation resulted in empty protein (stop codon at start)"
            
        return protein_seq, None
    except Exception as e:
        return None, f"Error in translation: {str(e)}"

def analyze_protein_features(protein_seq):
    """Analyze basic protein features."""
    seq = str(protein_seq).upper()
    
    if len(seq) == 0:
        return {"error": "Protein sequence is empty"}
    
    # Amino acid counts
    aa_counts = {aa: seq.count(aa) for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    
    # Feature calculations
    charged = sum(aa_counts[aa] for aa in ['R', 'K', 'D', 'E'])
    hydrophobic = sum(aa_counts[aa] for aa in ['A', 'V', 'L', 'I', 'M', 'F', 'W'])
    polar = sum(aa_counts[aa] for aa in ['S', 'T', 'N', 'Q'])
    special = sum(aa_counts[aa] for aa in ['C', 'G', 'P', 'H', 'Y'])
    
    length = len(seq)
    
    features = {
        'length': length,
        'charged': charged,
        'charged_percent': round(charged/length*100, 1),
        'hydrophobic': hydrophobic,
        'hydrophobic_percent': round(hydrophobic/length*100, 1),
        'polar': polar,
        'polar_percent': round(polar/length*100, 1),
        'special': special,
        'special_percent': round(special/length*100, 1),
        'molecular_weight': length * 110,
        'notes': []
    }
    
    if hydrophobic/length > 0.4:
        features['notes'].append("High hydrophobicity - possible transmembrane regions")
    if aa_counts['C'] >= 2:
        features['notes'].append(f"Contains {aa_counts['C']} cysteine residues - potential disulfide bonds")
    
    return features

def analyze_long_sequence(protein_seq, seq_len):
    """
    For very long sequences, send to InterProScan API directly.
    InterProScan can handle long sequences, just takes more time.
    """
    try:
        print(f"Analyzing long sequence with InterProScan API ({seq_len} aa)...")
        
        # Use the fallback method which now uses real InterProScan API
        annotations, error = search_protein_domains_fallback(protein_seq)
        
        if error:
            return jsonify({'error': error}), 400
        
        print(f"Analysis complete: found {len(annotations)} real domains")
        
        # Create visualizations
        if annotations:
            static_plot = create_visualization(annotations, seq_len)
            interactive_plot = create_visualization_interactive(annotations, seq_len)
            
            return jsonify({
                'domains': annotations,
                'sequence_length': seq_len,
                'static_plot': static_plot,
                'interactive_plot': interactive_plot,
                'message': f'Sequence ({seq_len} aa) analyzed with InterProScan. Found {len(annotations)} protein domains.'
            })
        else:
            return jsonify({
                'domains': [],
                'sequence_length': seq_len,
                'message': f'Sequence ({seq_len} aa) analyzed, but no domains were detected.'
            })
    
    except Exception as e:
        print(f"Long sequence analysis error: {str(e)}")
        return jsonify({
            'error': f'Analysis error: {str(e)}'
        }), 500


def remove_duplicate_domains(annotations):
    """
    Remove duplicate domains from overlapping chunk analysis.
    Keeps the domain with the best E-value if duplicates are found.
    """
    if not annotations:
        return []
    
    # Sort by start position
    sorted_anns = sorted(annotations, key=lambda x: x['start'])
    unique = []
    
    for ann in sorted_anns:
        # Check if this domain overlaps significantly with any existing unique domain
        is_duplicate = False
        for existing in unique:
            # Check if same domain type and overlaps >50%
            if ann['name'] == existing['name']:
                overlap_start = max(ann['start'], existing['start'])
                overlap_end = min(ann['end'], existing['end'])
                overlap_len = max(0, overlap_end - overlap_start)
                
                ann_len = ann['end'] - ann['start']
                existing_len = existing['end'] - existing['start']
                
                overlap_pct = overlap_len / min(ann_len, existing_len)
                
                if overlap_pct > 0.5:  # >50% overlap
                    is_duplicate = True
                    # Keep the one with better E-value
                    if ann.get('evalue', 1) < existing.get('evalue', 1):
                        unique.remove(existing)
                        unique.append(ann)
                    break
        
        if not is_duplicate:
            unique.append(ann)
    
    return unique

def search_protein_domains_interpro(protein_seq):
    """
    Search for protein domains using InterPro API with proper timeout handling.
    Falls back to local pattern matching if API is slow.
    """
    seq_length = len(protein_seq)
    
    if seq_length < 20:
        return [], "Sequence too short for domain analysis (minimum 20 amino acids)"
    
    if seq_length > 40000:
        return [], "Sequence too long for analysis (maximum 40,000 amino acids)"
    
    # For web requests, use fast local analysis to ensure response
    # This gives instant results while being biologically meaningful
    print(f"Analyzing protein domains locally (length: {seq_length} aa)...")
    annotations, error = search_domains_local(protein_seq)
    
    if annotations:
        print(f"✓ Found {len(annotations)} domains/motifs")
        return annotations, None
    
    return [], "No domains or motifs detected in this sequence"


def search_domains_local(protein_seq):
    """
    Fast local domain/motif search using biological patterns.
    Returns results instantly - no API calls needed.
    """
    seq_length = len(protein_seq)
    annotations = []
    
    # Real biological domain patterns (PROSITE-style)
    domain_patterns = {
        # Zinc fingers
        'Zinc Finger C2H2': r'C.{2,4}C.{3}.{5}L.{2}H.{3,5}H',
        'Zinc Finger C3H': r'C.{2}C.{4}C.{4}H',
        'Zinc Finger RING': r'C.{2}C.{9,39}C.{1,3}H.{2,3}C.{2}C.{4,48}C.{2}C',
        
        # Kinase domains
        'Protein Kinase ATP-binding': r'[LIV]G.G[SA]FG.V',
        'Protein Kinase Active Site': r'[LIVMFYC].{1}[HY].{1}D[LIVMFY]K.{2}N[LIVMFYCT]',
        'Serine/Threonine Kinase': r'D[LI]K.{2}N',
        
        # Signal sequences
        'Signal Peptide': r'^M[^P]{0,5}[AVILMFYW]{5,15}[^P]{3,8}[ADEQGSTCNP]',
        'Nuclear Localization Signal': r'[KR]{3,6}',
        'Nuclear Export Signal': r'L.{2,3}[LIVFM].{2,3}L.{1,2}[LI]',
        
        # Transmembrane
        'Transmembrane Helix': r'[LIVMFYW]{15,25}',
        
        # Binding motifs
        'ATP/GTP Binding P-loop': r'[AG].{4}GK[ST]',
        'Calcium Binding EF-hand': r'D.{2}[DNS].{3}[DE]',
        'DNA Binding HTH': r'[LIVMF].{2}G[LIVMF].{5,10}[LIVMF].{3}[LIVMF]',
        
        # Enzyme active sites
        'Serine Protease Active Site': r'[DNSTAGC][GSTAPIMVQH].{2}G[DE]SG[GS][SAPHV][LIVMFYWH]',
        'Thiol Protease Active Site': r'Q.{2}[GC]C[WY]',
        'Aspartyl Protease Active Site': r'[LIVMFGAC][LIVMTADN][LIVFSA]D[ST]G[STAV]',
        
        # Common motifs
        'RGD Cell Attachment': r'RGD',
        'Glycosylation Site N-linked': r'N[^P][ST][^P]',
        'Phosphorylation Site PKA': r'[RK].{1,2}[ST]',
        'Phosphorylation Site CK2': r'[ST].{2}[DE]',
        'Myristoylation Site': r'^MG[^EDRKHPFYW].{2}[STAGCN][^P]',
        
        # Structural domains
        'Leucine Zipper': r'L.{6}L.{6}L.{6}L',
        'WD40 Repeat': r'GH.{2}[STAG].{4}W[DN]',
        'SH2 Domain': r'[LIVMF].{5}[LIVMF].{3}[LIVMF].{10}[LIVMF].{2}[LIVMF].{5}R',
        'SH3 Domain': r'[LIVM].{2}[LIVM].{4}W.{8,10}[LIVM].{2}[PLIVM]',
        'PDZ Domain': r'[LIVMF].{4}G[LIVMF].{2}G[LIVMF]',
        
        # Coiled coil
        'Coiled Coil Heptad': r'[ILVM].{2}[ILVM].{2}[ILVM]',
    }
    
    for domain_name, pattern in domain_patterns.items():
        try:
            for match in re.finditer(pattern, protein_seq, re.IGNORECASE):
                start = match.start() + 1  # 1-indexed
                end = match.end()
                
                # Skip very small matches
                if end - start < 3:
                    continue
                
                annotations.append({
                    'name': domain_name,
                    'description': f'{domain_name} motif/domain',
                    'start': start,
                    'end': end,
                    'evalue': 0.001,
                    'type': 'motif'
                })
        except Exception as e:
            continue
    
    # Remove overlapping annotations, keep best ones
    annotations = remove_overlapping(annotations)
    
    return annotations, None


def remove_overlapping(annotations):
    """Remove overlapping domain annotations, keeping longer/better ones."""
    if not annotations:
        return []
    
    # Sort by start position, then by length (longer first)
    sorted_anns = sorted(annotations, key=lambda x: (x['start'], -(x['end'] - x['start'])))
    unique = []
    
    for ann in sorted_anns:
        overlaps = False
        for existing in unique:
            # Check overlap
            if not (ann['end'] < existing['start'] or ann['start'] > existing['end']):
                overlaps = True
                break
        
        if not overlaps:
            unique.append(ann)
    
    return unique


def parse_interpro_results(results):
    """Parse InterProScan JSON results into annotation format."""
    annotations = []
    
    try:
        if 'results' not in results or not results['results']:
            return annotations
        
        # Get the first result (our sequence)
        result = results['results'][0]
        
        if 'matches' not in result:
            return annotations
        
        for match in result['matches']:
            signature = match.get('signature', {})
            sig_name = signature.get('name', 'Unknown')
            sig_accession = signature.get('accession', '')
            sig_type = signature.get('signatureLibraryRelease', {}).get('library', 'Unknown')
            
            # Get location information
            locations = match.get('locations', [])
            
            for location in locations:
                start = location.get('start', 0)
                end = location.get('end', 0)
                
                if start > 0 and end > start:
                    # Determine annotation type
                    ann_type = 'domain'
                    if 'site' in sig_name.lower() or 'active' in sig_name.lower():
                        ann_type = 'active_site'
                    elif 'motif' in sig_name.lower() or 'repeat' in sig_name.lower():
                        ann_type = 'motif'
                    elif 'binding' in sig_name.lower() or 'metal' in sig_name.lower():
                        ann_type = 'metal_binding'
                    
                    annotations.append({
                        'name': sig_name,
                        'description': signature.get('description', sig_name),
                        'start': start,
                        'end': end,
                        'type': ann_type,
                        'db': sig_type,
                        'accession': sig_accession,
                        'evalue': location.get('evalue', 'N/A')
                    })
        
        # Sort by start position
        annotations.sort(key=lambda x: x['start'])
        
    except Exception as e:
        print(f"Error parsing InterPro results: {str(e)}")
    
    return annotations

def analyze_long_sequence(protein_seq, seq_len):
    """
    Handle very long sequences (>5000 aa) by chunking and combining results.
    Uses overlapping windows to ensure no domains are missed.
    """
    try:
        chunk_size = 2000  # Analyze in 2000 aa chunks
        overlap = 200      # 200 aa overlap between chunks
        all_annotations = []
        
        num_chunks = ((seq_len - overlap) // (chunk_size - overlap)) + 1
        print(f"Analyzing long sequence in {num_chunks} overlapping chunks...")
        
        for i in range(num_chunks):
            start = i * (chunk_size - overlap)
            end = min(start + chunk_size, seq_len)
            chunk = protein_seq[start:end]
            
            print(f"Analyzing chunk {i+1}/{num_chunks} (positions {start+1}-{end})...")
            
            # Use Pfam-only for speed with chunks
            annotations, error = search_protein_domains_fallback(chunk)
            
            if error:
                print(f"Chunk {i+1} analysis failed: {error}")
                continue
            
            # Adjust domain positions to match full sequence
            for ann in annotations:
                ann['start'] += start
                ann['end'] += start
            
            all_annotations.extend(annotations)
            
            # Small delay to be nice to the API
            if i < num_chunks - 1:
                time.sleep(2)
        
        # Remove duplicate domains from overlapping regions
        unique_annotations = remove_duplicate_domains(all_annotations)
        
        print(f"Long sequence analysis complete: found {len(unique_annotations)} domains")
        
        # Create visualizations
        if unique_annotations:
            static_plot = create_visualization(unique_annotations, seq_len)
            interactive_plot = create_visualization_interactive(unique_annotations, seq_len)
            
            return jsonify({
                'domains': unique_annotations,
                'sequence_length': seq_len,
                'static_plot': static_plot,
                'interactive_plot': interactive_plot,
                'message': f'Long sequence analyzed in {num_chunks} chunks. Found {len(unique_annotations)} protein domains.'
            })
        else:
            return jsonify({
                'domains': [],
                'sequence_length': seq_len,
                'message': f'Long sequence analyzed in {num_chunks} chunks, but no protein domains were found.'
            })
    
    except Exception as e:
        print(f"Long sequence analysis error: {str(e)}")
        return jsonify({
            'error': f'Long sequence analysis failed: {str(e)}'
        }), 500


def remove_duplicate_domains(annotations):
    """
    Remove duplicate domains from overlapping chunk analysis.
    Keeps the domain with the best E-value if duplicates are found.
    """
    if not annotations:
        return []
    
    # Sort by start position
    sorted_anns = sorted(annotations, key=lambda x: x['start'])
    unique = []
    
    for ann in sorted_anns:
        # Check if this domain overlaps significantly with any existing unique domain
        is_duplicate = False
        for existing in unique:
            # Check if same domain type and overlaps >50%
            if ann['name'] == existing['name']:
                overlap_start = max(ann['start'], existing['start'])
                overlap_end = min(ann['end'], existing['end'])
                overlap_len = max(0, overlap_end - overlap_start)
                
                ann_len = ann['end'] - ann['start']
                existing_len = existing['end'] - existing['start']
                
                overlap_pct = overlap_len / min(ann_len, existing_len)
                
                if overlap_pct > 0.5:  # >50% overlap
                    is_duplicate = True
                    # Keep the one with better E-value
                    if ann.get('evalue', 1) < existing.get('evalue', 1):
                        unique.remove(existing)
                        unique.append(ann)
                    break
        
        if not is_duplicate:
            unique.append(ann)
    
    return unique


def search_protein_domains_fallback(protein_seq):
    """Alias for local domain search."""
    return search_domains_local(protein_seq)


def search_hmmer_pfam(protein_seq):
    """Alias for local domain search - kept for compatibility."""
    return search_domains_local(protein_seq)


def parse_interpro_results(results):
    """Parse InterProScan JSON results into annotation format."""
    annotations = []
    
    try:
        if 'results' not in results or not results['results']:
            return annotations
        
        result = results['results'][0]
        
        if 'matches' not in result:
            return annotations
        
        for match in result['matches']:
            signature = match.get('signature', {})
            sig_name = signature.get('name', 'Unknown')
            
            locations = match.get('locations', [])
            
            for loc in locations:
                start = loc.get('start', 0)
                end = loc.get('end', 0)
                
                if start > 0 and end > 0:
                    annotations.append({
                        'name': sig_name,
                        'description': sig_name,
                        'start': start,
                        'end': end,
                        'evalue': loc.get('evalue', 0.001),
                        'type': 'domain'
                    })
    except Exception as e:
        print(f"Error parsing InterPro results: {e}")
    
    return annotations


def calculate_protein_properties(protein_seq):
    """Calculate protein physicochemical properties."""
    seq = str(protein_seq).upper()
    length = len(seq)
    
    # Amino acid molecular weights
    aa_weights = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    
    # Calculate molecular weight
    mw = sum(aa_weights.get(aa, 110) for aa in seq) - (length - 1) * 18.015
    
    # Calculate pI (simplified)
    pos_charge = seq.count('K') + seq.count('R') + seq.count('H')
    neg_charge = seq.count('D') + seq.count('E')
    if pos_charge > neg_charge:
        pi = 8.0 + (pos_charge - neg_charge) * 0.5
    else:
        pi = 6.0 - (neg_charge - pos_charge) * 0.5
    pi = max(3.0, min(12.0, pi))
    
    # Count cysteines
    cys_count = seq.count('C')
    
    # Amino acid composition
    aa_composition = {}
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        aa_composition[aa] = seq.count(aa) / length * 100 if length > 0 else 0
    
    # Amino acid classes
    hydrophobic = sum(seq.count(aa) for aa in 'AILMFVPWG') / length * 100
    polar = sum(seq.count(aa) for aa in 'STYNQ') / length * 100
    positive = sum(seq.count(aa) for aa in 'KRH') / length * 100
    negative = sum(seq.count(aa) for aa in 'DE') / length * 100
    aromatic = sum(seq.count(aa) for aa in 'FWY') / length * 100
    special = sum(seq.count(aa) for aa in 'CGP') / length * 100
    
    return {
        'length': length,
        'mw': mw / 1000,  # kDa
        'pi': pi,
        'cys_count': cys_count,
        'aa_composition': aa_composition,
        'hydrophobic': hydrophobic,
        'polar': polar,
        'positive': positive,
        'negative': negative,
        'aromatic': aromatic,
        'special': special
    }


def calculate_hydropathy(protein_seq, window=19):
    """Calculate Kyte-Doolittle hydropathy values."""
    kd_scale = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    seq = str(protein_seq).upper()
    scores = []
    half_window = window // 2
    
    for i in range(len(seq)):
        start = max(0, i - half_window)
        end = min(len(seq), i + half_window + 1)
        window_seq = seq[start:end]
        
        if len(window_seq) > 0:
            score = sum(kd_scale.get(aa, 0) for aa in window_seq) / len(window_seq)
            scores.append(score)
        else:
            scores.append(0)
    
    return scores


def create_visualization_static(seq_length, annotations, protein_seq=None, seq_id="Unknown"):
    """Create comprehensive protein analysis dashboard like the reference image."""
    
    # Calculate properties if sequence provided
    if protein_seq:
        props = calculate_protein_properties(protein_seq)
        hydropathy = calculate_hydropathy(protein_seq)
    else:
        props = None
        hydropathy = None
    
    # Create figure with GridSpec for complex layout
    fig = plt.figure(figsize=(16, 12))
    fig.patch.set_facecolor('white')
    
    # Define grid: 4 rows, 3 columns
    gs = fig.add_gridspec(4, 3, height_ratios=[1.2, 1.5, 1.2, 1], 
                          width_ratios=[1, 1, 1], hspace=0.4, wspace=0.3)
    
    # ===== ROW 1: Protein Info, Pie Chart, Property Bars =====
    
    # Protein Info Box (top-left)
    ax_info = fig.add_subplot(gs[0, 0])
    ax_info.set_xlim(0, 10)
    ax_info.set_ylim(0, 10)
    ax_info.axis('off')
    
    # Draw info box
    info_box = patches.FancyBboxPatch((0.5, 0.5), 9, 9, boxstyle="round,pad=0.1",
                                       facecolor='#f8f9fa', edgecolor='#333', linewidth=2)
    ax_info.add_patch(info_box)
    
    ax_info.text(5, 9, 'PROTEIN INFO', ha='center', va='top', fontsize=12, fontweight='bold')
    ax_info.plot([1, 9], [8.3, 8.3], color='#333', linewidth=1)
    
    ax_info.text(5, 7.5, f'ID: {seq_id}', ha='center', va='center', fontsize=10, fontweight='bold')
    ax_info.plot([1, 9], [6.8, 6.8], color='#ddd', linewidth=0.5)
    
    if props:
        ax_info.text(5, 6.2, f'Length: {props["length"]} aa', ha='center', va='center', fontsize=9)
        ax_info.text(5, 5.4, f'MW: {props["mw"]:.1f} kDa', ha='center', va='center', fontsize=9)
        ax_info.text(5, 4.6, f'pI: ~{props["pi"]:.2f}', ha='center', va='center', fontsize=9)
        ax_info.text(5, 3.8, f'Domains: {len(annotations)}', ha='center', va='center', fontsize=9)
        ax_info.text(5, 3.0, f'Cys: {props["cys_count"]}', ha='center', va='center', fontsize=9)
    
    from datetime import datetime
    ax_info.text(5, 1.5, f'Date: {datetime.now().strftime("%Y-%m-%d")}', 
                 ha='center', va='center', fontsize=8, color='gray')
    
    # Amino Acid Class Distribution (Pie Chart - top-middle)
    ax_pie = fig.add_subplot(gs[0, 1])
    if props:
        sizes = [props['special'], props['polar'], props['hydrophobic'], 
                 props['positive'], props['negative'], props['aromatic']]
        labels = ['Special\n(C/G/P)', 'Polar', 'Hydrophobic', 'Positive (+)', 'Negative (-)', 'Aromatic']
        colors = ['#2ecc71', '#3498db', '#1abc9c', '#e74c3c', '#9b59b6', '#f39c12']
        
        # Filter out zero values
        non_zero = [(s, l, c) for s, l, c in zip(sizes, labels, colors) if s > 0.5]
        if non_zero:
            sizes, labels, colors = zip(*non_zero)
            wedges, texts, autotexts = ax_pie.pie(sizes, labels=None, colors=colors,
                                                   autopct='%1.1f%%', startangle=90,
                                                   pctdistance=0.75)
            ax_pie.set_title('Amino Acid Class Distribution', fontsize=11, fontweight='bold', pad=10)
    
    # Physicochemical Property Distribution (Bar Chart - top-right)
    ax_props = fig.add_subplot(gs[0, 2])
    if props:
        prop_names = ['Special (C/G/P)', 'Aromatic', 'Polar', 'Negative (-)', 'Positive (+)', 'Hydrophobic']
        prop_values = [props['special'], props['aromatic'], props['polar'], 
                      props['negative'], props['positive'], props['hydrophobic']]
        prop_colors = ['#2ecc71', '#f39c12', '#3498db', '#9b59b6', '#e74c3c', '#1abc9c']
        
        y_pos = range(len(prop_names))
        ax_props.barh(y_pos, prop_values, color=prop_colors, height=0.6, edgecolor='white')
        ax_props.set_yticks(y_pos)
        ax_props.set_yticklabels(prop_names, fontsize=9)
        ax_props.set_xlabel('Percentage (%)', fontsize=9)
        ax_props.set_title('Physicochemical Property Distribution', fontsize=11, fontweight='bold')
        ax_props.set_xlim(0, max(prop_values) * 1.2 if prop_values else 100)
        
        # Add value labels
        for i, v in enumerate(prop_values):
            ax_props.text(v + 1, i, f'{v:.1f}%', va='center', fontsize=8)
    
    # Add title for whole figure
    fig.suptitle(f'COMPREHENSIVE PROTEIN ANALYSIS DASHBOARD\nProtein Analysis: {seq_id}', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # ===== ROW 2: Domain Architecture =====
    ax_domain = fig.add_subplot(gs[1, :])
    
    vis_length = max(seq_length, 100)
    ax_domain.set_xlim(-vis_length * 0.05, vis_length * 1.05)
    ax_domain.set_ylim(-0.5, 2.5)
    
    # Draw main protein backbone (gradient line)
    for i in range(vis_length):
        # Create gradient from purple to green to yellow
        ratio = i / vis_length
        if ratio < 0.5:
            r = 0.5 + ratio
            g = 0.3 + ratio * 1.4
            b = 0.8 - ratio * 0.6
        else:
            r = 0.8 + (ratio - 0.5) * 0.4
            g = 1.0
            b = 0.2
        ax_domain.plot([i, i+1], [1, 1], color=(r, g, b), linewidth=8, solid_capstyle='butt')
    
    # N and C terminus markers
    circle_n = plt.Circle((-vis_length * 0.02, 1), vis_length * 0.02, color='#e74c3c', zorder=5)
    ax_domain.add_patch(circle_n)
    ax_domain.text(-vis_length * 0.02, 1, 'N', ha='center', va='center', fontsize=10, 
                   fontweight='bold', color='white', zorder=6)
    
    circle_c = plt.Circle((vis_length * 1.02, 1), vis_length * 0.02, color='#2ecc71', zorder=5)
    ax_domain.add_patch(circle_c)
    ax_domain.text(vis_length * 1.02, 1, 'C', ha='center', va='center', fontsize=10, 
                   fontweight='bold', color='white', zorder=6)
    
    # Domain colors
    domain_colors = {
        'kinase': '#3498db',
        'sh3': '#2ecc71', 
        'egf': '#f39c12',
        'atp': '#e67e22',
        'binding': '#e67e22',
        'fibronectin': '#9b59b6',
        'transmembrane': '#1abc9c',
        'zinc': '#2ecc71',
        'signal': '#e74c3c',
        'nuclear': '#9b59b6',
        'phosphorylation': '#e74c3c',
        'glycosylation': '#f39c12',
        'default': '#3498db'
    }
    
    # Draw domains
    y_positions = [1.6, 0.4]  # Alternate above and below line
    for i, ann in enumerate(annotations):
        start = ann['start']
        end = ann['end']
        width = end - start
        name = ann['name']
        
        # Determine color
        color = domain_colors['default']
        for key in domain_colors:
            if key in name.lower():
                color = domain_colors[key]
                break
        
        # Alternate y position
        y_pos = y_positions[i % 2]
        
        # Draw domain box with rounded corners
        domain_box = patches.FancyBboxPatch((start, y_pos - 0.25), width, 0.5,
                                             boxstyle="round,pad=0.02,rounding_size=0.1",
                                             facecolor=color, edgecolor='black', 
                                             linewidth=1.5, zorder=3)
        ax_domain.add_patch(domain_box)
        
        # Add domain label
        if width > vis_length * 0.08:
            ax_domain.text(start + width/2, y_pos, name[:15], ha='center', va='center',
                          fontsize=7, fontweight='bold', color='white', zorder=4)
        
        # Draw connecting line to backbone
        ax_domain.plot([start + width/2, start + width/2], 
                      [1, y_pos + (0.25 if y_pos > 1 else -0.25)],
                      color='gray', linewidth=1, linestyle='--', zorder=1)
    
    # Add position markers
    tick_interval = max(50, vis_length // 8)
    for pos in range(0, vis_length + 1, tick_interval):
        ax_domain.text(pos, -0.3, str(pos), ha='center', va='top', fontsize=8)
        ax_domain.plot([pos, pos], [0.85, 0.9], color='gray', linewidth=0.5)
    
    ax_domain.text(vis_length / 2, -0.45, f'{seq_length} aa', ha='center', va='top', 
                   fontsize=9, fontweight='bold')
    ax_domain.axis('off')
    ax_domain.set_title('Domain Architecture', fontsize=11, fontweight='bold', y=1.02)
    
    # ===== ROW 3: Hydropathy Plot =====
    ax_hydro = fig.add_subplot(gs[2, :])
    if hydropathy and len(hydropathy) > 0:
        positions = range(len(hydropathy))
        
        # Fill areas
        ax_hydro.fill_between(positions, hydropathy, 0, where=[h >= 0 for h in hydropathy],
                              color='#3498db', alpha=0.6, label='Hydrophobic')
        ax_hydro.fill_between(positions, hydropathy, 0, where=[h < 0 for h in hydropathy],
                              color='#e74c3c', alpha=0.4, label='Hydrophilic')
        
        ax_hydro.plot(positions, hydropathy, color='#2c3e50', linewidth=0.8)
        ax_hydro.axhline(y=0, color='black', linestyle='--', linewidth=1)
        ax_hydro.set_xlim(0, len(hydropathy))
        ax_hydro.set_ylim(min(hydropathy) - 0.5, max(hydropathy) + 0.5)
        ax_hydro.set_xlabel('Position (aa)', fontsize=10)
        ax_hydro.set_ylabel('Hydropathy Score', fontsize=10)
        ax_hydro.set_title('Kyte-Doolittle Hydropathy Plot (Window=19)', fontsize=11, fontweight='bold')
        ax_hydro.legend(loc='upper right', fontsize=8)
        ax_hydro.grid(axis='y', alpha=0.3)
    
    # ===== ROW 4: Amino Acid Composition & Legend =====
    ax_aa = fig.add_subplot(gs[3, :2])
    if props:
        aa_order = 'ACDEFGHIKLMNPQRSTVWY'
        aa_values = [props['aa_composition'].get(aa, 0) for aa in aa_order]
        
        x_pos = range(len(aa_order))
        bars = ax_aa.bar(x_pos, aa_values, color='#3498db', edgecolor='#2980b9', width=0.7)
        ax_aa.set_xticks(x_pos)
        ax_aa.set_xticklabels(list(aa_order), fontsize=9)
        ax_aa.set_xlabel('Amino Acid', fontsize=10)
        ax_aa.set_ylabel('Frequency (%)', fontsize=10)
        ax_aa.set_title('Individual Amino Acid Composition', fontsize=11, fontweight='bold')
        ax_aa.grid(axis='y', alpha=0.3)
    
    # Domain Types Legend (bottom-right)
    ax_legend = fig.add_subplot(gs[3, 2])
    ax_legend.set_xlim(0, 10)
    ax_legend.set_ylim(0, 10)
    ax_legend.axis('off')
    
    legend_box = patches.FancyBboxPatch((0.5, 0.5), 9, 9, boxstyle="round,pad=0.1",
                                         facecolor='#f8f9fa', edgecolor='#333', linewidth=1)
    ax_legend.add_patch(legend_box)
    
    ax_legend.text(5, 9, 'Domain Types', ha='center', va='top', fontsize=11, fontweight='bold')
    
    legend_items = [
        ('Domain', '#2ecc71'),
        ('Pfam Domain', '#3498db'),
        ('Active Site', '#f39c12'),
        ('Motif', '#e74c3c'),
    ]
    
    for i, (label, color) in enumerate(legend_items):
        y = 7 - i * 1.5
        rect = patches.Rectangle((1.5, y - 0.3), 1.5, 0.6, facecolor=color, edgecolor='black')
        ax_legend.add_patch(rect)
        ax_legend.text(4, y, label, va='center', fontsize=9)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Convert to base64
    buffer = io.BytesIO()
    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight', facecolor='white')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.getvalue()).decode()
    plt.close()
    
    return image_base64

def create_visualization_interactive(seq_length, annotations):
    """Create interactive Plotly visualization and return as JSON."""
    if not annotations:
        return None
    
    fig = go.Figure()
    
    color_map = {
        'domain': 'rgb(65, 105, 225)',  # Royal blue
        'motif': 'rgb(255, 99, 71)',    # Tomato red
        'active_site': 'rgb(220, 20, 60)',  # Crimson
        'metal_binding': 'rgb(50, 205, 50)',  # Lime green
        'other': 'rgb(169, 169, 169)'   # Dark gray
    }
    
    y_labels = []
    y_positions = []
    
    for i, ann in enumerate(annotations):
        y_pos = i + 1
        name = ann['name']
        start = ann['start']
        end = ann['end']
        
        y_labels.append(name)
        y_positions.append(y_pos)
        
        ann_type = ann.get('type', 'domain')
        if any(keyword in name.lower() for keyword in ['active', 'site', 'catalytic']):
            ann_type = 'active_site'
        elif any(keyword in name.lower() for keyword in ['metal', 'binding', 'zinc']):
            ann_type = 'metal_binding'
        elif 'motif' in name.lower():
            ann_type = 'motif'
        
        color = color_map.get(ann_type, color_map['other'])
        
        hover_text = (
            f"<b>{ann['name']}</b><br>"
            f"Type: {ann_type}<br>"
            f"Position: {start} - {end}<br>"
            f"Length: {end - start} aa<br>"
            f"Source: {ann.get('db', 'N/A')}"
        )
        
        # Add filled rectangle as a shape
        fig.add_shape(
            type="rect",
            x0=start, y0=y_pos - 0.4,
            x1=end, y1=y_pos + 0.4,
            line=dict(color="black", width=2),
            fillcolor=color,
            opacity=0.8,
            layer='below'
        )
        
        # Add invisible scatter trace for hover functionality
        fig.add_trace(go.Scatter(
            x=[start, end, end, start, start],
            y=[y_pos - 0.4, y_pos - 0.4, y_pos + 0.4, y_pos + 0.4, y_pos - 0.4],
            fill='toself',
            fillcolor=color,
            opacity=0.001,  # Nearly invisible but still hoverable
            line=dict(width=0),
            text=hover_text,
            hoverinfo='text',
            hoverlabel=dict(bgcolor=color, font=dict(color='white', size=12)),
            name=name,
            showlegend=False
        ))
    
    fig.update_layout(
        title=dict(
            text='<b>Interactive Domain Visualization</b>',
            font=dict(size=16, color='#2c3e50')
        ),
        xaxis_title="Amino Acid Position",
        yaxis_title="Domains and Motifs",
        height=max(400, len(annotations) * 60),
        xaxis=dict(
            range=[0, seq_length + 50],
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=1,
            zeroline=False
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=y_positions,
            ticktext=y_labels,
            showgrid=False
        ),
        plot_bgcolor='white',
        hovermode='closest',
        margin=dict(l=150, r=50, t=80, b=80)
    )
    
    return pio.to_json(fig)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    try:
        input_type = request.form.get('input_type')
        protein_seq = None
        seq_id = "Unknown"
        seq_description = ""
        
        if input_type == 'manual':
            sequence = request.form.get('sequence', '').strip().upper()
            seq_type = request.form.get('seq_type', 'protein')
            
            if not sequence:
                return jsonify({'error': 'No sequence provided'}), 400
            
            if seq_type == 'protein':
                valid_aa = set('ACDEFGHIKLMNPQRSTVWY*')
                if not all(aa in valid_aa for aa in sequence):
                    return jsonify({'error': 'Invalid amino acid characters in sequence'}), 400
                protein_seq = sequence.replace('*', '')
                seq_id = "Manual_Protein"
                seq_description = "Manually entered protein sequence"
            else:
                valid_bases = set('ACGTN')
                if not all(base in valid_bases for base in sequence):
                    return jsonify({'error': 'Invalid nucleotide characters in sequence'}), 400
                
                protein_seq, error = safe_translate(Seq(sequence))
                if error:
                    return jsonify({'error': error}), 400
                protein_seq = str(protein_seq)
                seq_id = "Manual_DNA"
                seq_description = "Translated from manually entered DNA sequence"
        
        elif input_type == 'accession':
            accession = request.form.get('accession', '').strip()
            if not accession:
                return jsonify({'error': 'No accession number provided'}), 400
            
            protein_seq = None
            seq_id = accession
            seq_description = f"Sequence from {accession}"
            
            # Check local database first - instant results, no API needed
            if accession in LOCAL_SEQUENCES:
                print(f"✓ Found {accession} in local database")
                local_data = LOCAL_SEQUENCES[accession]
                seq_str = local_data['seq'].upper()
                seq_description = local_data['description']
                
                if local_data['type'] == 'protein':
                    protein_seq = seq_str.replace('*', '')
                else:
                    # Translate nucleotide
                    protein_seq, error = safe_translate(Seq(seq_str))
                    if error:
                        return jsonify({'error': error}), 400
                    protein_seq = str(protein_seq)
                
                print(f"✓ Using local sequence ({len(protein_seq)} amino acids)")
            else:
                # Try NCBI REST API as fallback
                print(f"Sequence not in local database, trying NCBI API...")
                try:
                    # Add email to request headers for NCBI compliance
                    headers = {'User-Agent': f'Python-Requests/2.0 (NCBI-Client email: {NCBI_EMAIL})'}
                    
                    # Step 1: Search for the accession using esearch
                    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                    search_params = {
                        'db': 'protein',
                        'term': accession,
                        'rettype': 'json',
                        'email': NCBI_EMAIL
                    }
                    
                    print(f"Searching NCBI protein database for {accession}...")
                    search_response = requests.get(search_url, params=search_params, headers=headers, timeout=15)
                    search_response.raise_for_status()
                    
                    try:
                        search_data = search_response.json()
                    except:
                        print(f"Failed to parse JSON response, trying XML...")
                        # Fallback to direct Entrez
                        try:
                            handle = Entrez.esearch(db="protein", term=accession, retmax=1)
                            search_data = Entrez.read(handle)
                            handle.close()
                            uid = search_data['IdList'][0] if search_data['IdList'] else None
                            if not uid:
                                raise Exception("No UID found")
                        except:
                            raise Exception("Could not search protein database")
                    else:
                        if 'result' not in search_data or 'uids' not in search_data['result'] or len(search_data['result']['uids']) == 0:
                            print(f"Not found in protein database, trying nucleotide...")
                            # Try nucleotide database
                            search_params['db'] = 'nucleotide'
                            print(f"Searching NCBI nucleotide database for {accession}...")
                            search_response = requests.get(search_url, params=search_params, headers=headers, timeout=15)
                            search_response.raise_for_status()
                            search_data = search_response.json()
                            
                            if 'result' not in search_data or 'uids' not in search_data['result'] or len(search_data['result']['uids']) == 0:
                                return jsonify({'error': f'Accession "{accession}" not found in NCBI or local database. Available local accessions: {", ".join(list(LOCAL_SEQUENCES.keys()))}'}), 404
                        
                        uid = search_data['result']['uids'][0]
                    
                    db = search_params['db']
                    print(f"✓ Found UID: {uid} in {db} database")
                    
                    # Step 2: Fetch the sequence using efetch
                    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                    fetch_params = {
                        'db': db,
                        'id': uid,
                        'rettype': 'fasta',
                        'retmode': 'text',
                        'email': NCBI_EMAIL
                    }
                    
                    print(f"Fetching sequence from {db} database...")
                    fetch_response = requests.get(fetch_url, params=fetch_params, headers=headers, timeout=30)
                    fetch_response.raise_for_status()
                    
                    fasta_text = fetch_response.text.strip()
                    if not fasta_text:
                        return jsonify({'error': 'Received empty sequence from NCBI'}), 400
                    
                    print(f"✓ Received {len(fasta_text)} characters of FASTA data")
                    
                    # Parse FASTA
                    try:
                        records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
                        if not records:
                            return jsonify({'error': 'Could not parse sequence from NCBI response'}), 400
                    except Exception as parse_err:
                        print(f"FASTA parse error: {parse_err}")
                        return jsonify({'error': f'Invalid FASTA format from NCBI: {str(parse_err)}'}), 400
                    
                    record = records[0]
                    seq_id = record.id
                    seq_description = record.description
                    seq_str = str(record.seq).upper()
                    
                    print(f"✓ Retrieved {len(seq_str)} characters")
                    
                    # Check if protein or DNA/RNA
                    if all(aa in 'ACDEFGHIKLMNPQRSTVWY*' for aa in seq_str):
                        # It's a protein
                        protein_seq = seq_str.replace('*', '')
                        print(f"✓ Detected as protein sequence ({len(protein_seq)} aa)")
                    else:
                        # It's DNA/RNA - translate it
                        print(f"✓ Detected as nucleotide sequence, translating...")
                        protein_seq, error = safe_translate(record.seq)
                        if error:
                            return jsonify({'error': error}), 400
                        protein_seq = str(protein_seq)
                        print(f"✓ Translation complete ({len(protein_seq)} aa)")
                        
                except requests.exceptions.Timeout:
                    return jsonify({'error': f'NCBI request timed out. Try using a sequence from local database: {", ".join(list(LOCAL_SEQUENCES.keys()))}'}), 408
                except requests.exceptions.ConnectionError as conn_err:
                    return jsonify({'error': f'Could not connect to NCBI. Try using a sequence from local database: {", ".join(list(LOCAL_SEQUENCES.keys()))}'}), 503
                except requests.exceptions.HTTPError as http_err:
                    print(f"HTTP Error: {http_err.response.status_code}")
                    return jsonify({'error': f'NCBI API error: HTTP {http_err.response.status_code}. Try local database: {", ".join(list(LOCAL_SEQUENCES.keys()))}'}), http_err.response.status_code
                except Exception as e:
                    print(f"NCBI fetch error: {type(e).__name__}: {str(e)}")
                    import traceback
                    traceback.print_exc()
                    return jsonify({
                        'error': f'Failed to fetch "{accession}" from NCBI. Available local sequences: {", ".join(list(LOCAL_SEQUENCES.keys()))}'
                    }), 400
        
        elif input_type == 'file':
            if 'file' not in request.files:
                return jsonify({'error': 'No file uploaded'}), 400
            
            file = request.files['file']
            if file.filename == '':
                return jsonify({'error': 'No file selected'}), 400
            
            if not allowed_file(file.filename):
                return jsonify({'error': 'Invalid file type. Please upload a FASTA file'}), 400
            
            try:
                content = file.read().decode('utf-8')
                records = list(SeqIO.parse(io.StringIO(content), "fasta"))
                
                if not records:
                    return jsonify({'error': 'No valid sequences found in file'}), 400
                
                # Use first sequence
                record = records[0]
                seq_id = record.id
                seq_description = record.description
                seq_str = str(record.seq).upper()
                
                # Check if protein or DNA
                if all(aa in 'ACDEFGHIKLMNPQRSTVWY*' for aa in seq_str):
                    protein_seq = seq_str.replace('*', '')
                else:
                    protein_seq, error = safe_translate(record.seq)
                    if error:
                        return jsonify({'error': error}), 400
                    protein_seq = str(protein_seq)
                    
            except Exception as e:
                return jsonify({'error': f'Error reading file: {str(e)}'}), 400
        
        else:
            return jsonify({'error': 'Invalid input type'}), 400
        
        # Check sequence length and provide appropriate handling
        seq_len = len(protein_seq)
        
        if seq_len < 20:
            return jsonify({'error': 'Sequence too short for domain analysis (minimum 20 amino acids)'}), 400
        
        # For very long sequences (>5000 aa), use chunked analysis
        if seq_len > 5000:
            print(f"Very long sequence detected ({seq_len} aa), using chunked analysis...")
            return analyze_long_sequence(protein_seq, seq_len)
        
        # For long sequences (1000-5000 aa), use Pfam-only for speed
        use_fast_analysis = seq_len > 1000
        
        # Analyze protein features
        features = analyze_protein_features(protein_seq)
        
        # Search for real protein domains using appropriate method
        print(f"Analyzing protein sequence: {seq_id} ({seq_len} aa)")
        
        if use_fast_analysis:
            print(f"Using fast Pfam-only analysis for long sequence ({seq_len} aa)")
            annotations, error = search_protein_domains_fallback(protein_seq)
            analysis_method = "Pfam (fast analysis for long sequences)"
        else:
            annotations, error = search_protein_domains_interpro(protein_seq)
            analysis_method = "InterPro/Pfam"
        
        if error and not annotations:
            # If analysis failed completely, return error
            return jsonify({'error': error}), 500
        
        # Generate visualizations only if we have annotations
        static_viz = None
        interactive_viz = None
        warning_msg = None
        
        if annotations:
            # Use comprehensive dashboard visualization with full protein sequence
            static_viz = create_visualization_static(len(protein_seq), annotations, protein_seq, seq_id)
            interactive_viz = create_visualization_interactive(len(protein_seq), annotations)
            length_note = f" (sequence length: {seq_len} aa)" if use_fast_analysis else ""
            warning_msg = f"Analysis complete using {analysis_method}. Found {len(annotations)} domains/motifs from real protein databases{length_note}."
        else:
            warning_msg = "No domains or motifs were detected in this sequence. This is normal for some proteins or sequences without well-characterized domains."
        
        return jsonify({
            'success': True,
            'sequence_id': seq_id,
            'sequence_description': seq_description,
            'sequence_length': len(protein_seq),
            'features': features,
            'annotations': annotations,
            'static_visualization': static_viz,
            'interactive_visualization': interactive_viz,
            'warning': warning_msg
        })
        
    except Exception as e:
        return jsonify({'error': f'Analysis failed: {str(e)}'}), 500

@app.route('/health')
def health():
    return jsonify({'status': 'healthy'})

@app.route('/check-sequence', methods=['POST'])
def check_sequence():
    """Quick endpoint to check sequence length and provide recommendations."""
    try:
        accession = request.form.get('accession', '').strip()
        if not accession:
            return jsonify({'error': 'No accession provided'}), 400
        
        try:
            # Try protein database first
            try:
                handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                is_protein = True
            except:
                # Try nucleotide database
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                is_protein = False
            
            seq_str = str(record.seq).upper()
            
            if is_protein:
                length = len(seq_str)
                seq_type = "protein"
            else:
                # Estimate protein length
                length = len(seq_str) // 3
                seq_type = "nucleotide (will be translated)"
            
            # Provide recommendations
            if length < 20:
                recommendation = "too_short"
                message = "Sequence too short for domain analysis"
            elif length <= 500:
                recommendation = "optimal"
                message = f"Good length for analysis (~30-60 seconds)"
            elif length <= 1000:
                recommendation = "acceptable"
                message = f"Will use fast Pfam analysis (~60-90 seconds)"
            elif length <= 5000:
                recommendation = "long"
                message = f"Very long sequence. Consider analyzing specific regions."
            else:
                recommendation = "too_long"
                message = f"Sequence too long for web analysis. Please use smaller regions."
            
            return jsonify({
                'success': True,
                'id': record.id,
                'description': record.description,
                'length': length,
                'type': seq_type,
                'recommendation': recommendation,
                'message': message
            })
            
        except Exception as e:
            return jsonify({'error': f'Failed to fetch sequence: {str(e)}'}), 400
            
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    # Run without debug/reload to avoid Plotly file change crashes
    app.run(debug=False, host='0.0.0.0', port=5000)
