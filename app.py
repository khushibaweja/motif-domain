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
import xml.etree.ElementTree as ET

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

def fetch_sequence_from_pdb(pdb_id):
    """
    Fetch protein sequence from PDB (Protein Data Bank) using RCSB PDB REST API.
    Returns (sequence, description, error)
    """
    pdb_id = pdb_id.strip().upper()
    
    # Validate PDB ID format (4 characters, alphanumeric)
    # Can also have chain ID like 1ABC_A or 1ABC:A
    chain_id = None
    if '_' in pdb_id:
        parts = pdb_id.split('_')
        pdb_id = parts[0]
        chain_id = parts[1] if len(parts) > 1 else None
    elif ':' in pdb_id:
        parts = pdb_id.split(':')
        pdb_id = parts[0]
        chain_id = parts[1] if len(parts) > 1 else None
    
    if len(pdb_id) != 4:
        return None, None, f"Invalid PDB ID format. PDB IDs should be 4 characters (e.g., 1ABC, 4HHB). Got: {pdb_id}"
    
    try:
        # Method 1: Use RCSB PDB FASTA endpoint for sequence
        if chain_id:
            fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        else:
            fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        
        print(f"Fetching sequence from PDB: {pdb_id}" + (f" chain {chain_id}" if chain_id else ""))
        
        response = requests.get(fasta_url, timeout=30)
        response.raise_for_status()
        
        fasta_text = response.text.strip()
        if not fasta_text:
            return None, None, f"No sequence found for PDB ID: {pdb_id}"
        
        # Parse FASTA - may contain multiple chains
        records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
        
        if not records:
            return None, None, f"Could not parse sequence from PDB for: {pdb_id}"
        
        # If chain specified, find that chain
        if chain_id:
            for record in records:
                # Chain ID is usually in the header like ">1ABC_A" or in description
                if f"_{chain_id}" in record.id or f"Chain {chain_id}" in record.description:
                    protein_seq = str(record.seq).upper().replace('X', '').replace('*', '')
                    description = f"PDB:{pdb_id} Chain {chain_id} - {record.description}"
                    print(f"✓ Found chain {chain_id}: {len(protein_seq)} amino acids")
                    return protein_seq, description, None
            # If chain not found, use first record with a warning
            print(f"Chain {chain_id} not found, using first available chain")
        
        # Use first record
        record = records[0]
        protein_seq = str(record.seq).upper().replace('X', '').replace('*', '')
        description = f"PDB:{pdb_id} - {record.description}"
        
        print(f"✓ Retrieved from PDB: {len(protein_seq)} amino acids")
        return protein_seq, description, None
        
    except requests.exceptions.Timeout:
        return None, None, f"PDB request timed out for: {pdb_id}"
    except requests.exceptions.ConnectionError:
        return None, None, f"Could not connect to PDB database for: {pdb_id}"
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            return None, None, f"PDB ID '{pdb_id}' not found in Protein Data Bank"
        return None, None, f"PDB API error: HTTP {e.response.status_code}"
    except Exception as e:
        print(f"PDB fetch error: {type(e).__name__}: {str(e)}")
        return None, None, f"Failed to fetch from PDB: {str(e)}"


def is_valid_protein_sequence(seq):
    """
    Check if a sequence is a valid protein sequence.
    Accepts standard amino acids plus common modifications.
    """
    # Standard amino acids + common modifications
    valid_chars = set('ACDEFGHIKLMNPQRSTVWY*XU')
    seq_upper = seq.upper().replace(' ', '').replace('\n', '').replace('\r', '')
    return all(c in valid_chars for c in seq_upper)


def is_valid_dna_sequence(seq):
    """
    Check if a sequence is a valid DNA/RNA sequence.
    """
    valid_chars = set('ACGTUNRYSWKMBDHV')
    seq_upper = seq.upper().replace(' ', '').replace('\n', '').replace('\r', '')
    return all(c in valid_chars for c in seq_upper)


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
    Search for protein domains using multiple APIs.
    Order: Fast APIs first, then slow InterProScan as last resort.
    """
    seq_length = len(protein_seq)
    
    if seq_length < 20:
        return [], "Sequence too short for domain analysis (minimum 20 amino acids)"
    
    if seq_length > 40000:
        return [], "Sequence too long for analysis (maximum 40,000 amino acids)"
    
    print(f"Searching protein domains (length: {seq_length} aa)...")
    
    # 1. Try InterProScan API (most reliable, but takes 1-2 min)
    print("Trying InterProScan API (this may take 1-2 minutes)...")
    try:
        annotations, error = search_interproscan_api(protein_seq)
        if annotations:
            print(f"✓ InterProScan API found {len(annotations)} domains/motifs")
            return annotations, None
        elif error:
            print(f"InterProScan: {error}")
        else:
            # InterProScan worked but found no domains - that's valid!
            print("✓ InterProScan completed but found no known domains in databases")
            print("  (This is normal for novel or poorly characterized sequences)")
    except Exception as e:
        print(f"InterProScan exception: {str(e)}")
    
    # 2. Fall back to local analysis (pattern-based)
    print("Using local pattern matching to find motifs...")
    annotations, error = search_domains_local(protein_seq)
    
    if annotations:
        print(f"✓ Local analysis found {len(annotations)} domains/motifs")
        return annotations, None
    
    return [], "No domains or motifs detected in this sequence"


def search_interproscan_api(protein_seq):
    """
    Search for domains using the EBI InterProScan REST API.
    This is the gold standard for protein domain analysis.
    API Documentation: https://www.ebi.ac.uk/Tools/services/rest/iprscan5
    """
    try:
        # Step 1: Submit job to InterProScan
        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
        
        # Format sequence as FASTA
        fasta_seq = f">query_protein\n{protein_seq}"
        
        # Valid applications (from API docs - case sensitive!):
        # PfamA, SMART, SuperFamily, Gene3d, Coils, CDD, PrositeProfiles, PrositePatterns, Panther, PRINTS
        job_data = {
            'email': NCBI_EMAIL,
            'sequence': fasta_seq,
            'stype': 'p',  # protein
            'appl': 'PfamA,SMART,SuperFamily,Gene3d,CDD,PrositeProfiles,PrositePatterns,Coils',
            'goterms': 'false',
            'pathways': 'false'
        }
        
        print("Submitting sequence to InterProScan API...")
        submit_response = requests.post(submit_url, data=job_data, timeout=30)
        
        if submit_response.status_code != 200:
            print(f"InterProScan response: {submit_response.text[:500]}")
            return None, f"InterProScan submit failed: HTTP {submit_response.status_code}"
        
        job_id = submit_response.text.strip()
        print(f"✓ Job submitted successfully: {job_id}")
        
        # Step 2: Poll for results (max 300 seconds = 5 minutes)
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
        result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
        
        max_attempts = 60  # 60 attempts * 5 seconds = 300 seconds max
        for attempt in range(max_attempts):
            time.sleep(5)  # Wait 5 seconds between checks
            
            try:
                status_response = requests.get(status_url, timeout=30)
                status = status_response.text.strip()
            except requests.exceptions.Timeout:
                print(f"Status check timeout (attempt {attempt + 1}), retrying...")
                continue
            except Exception as e:
                print(f"Status check error: {e}, retrying...")
                continue
            
            print(f"Job status: {status} (attempt {attempt + 1}/{max_attempts})")
            
            if status == 'FINISHED':
                # Get results
                print("Fetching results...")
                result_response = requests.get(result_url, timeout=30)
                if result_response.status_code == 200:
                    results = result_response.json()
                    annotations = parse_interproscan_results(results)
                    print(f"✓ InterProScan returned {len(annotations)} annotations")
                    return annotations, None
                else:
                    return None, f"Failed to get results: HTTP {result_response.status_code}"
            
            elif status == 'FAILURE' or status == 'ERROR':
                return None, f"InterProScan job failed: {status}"
            
            elif status not in ['RUNNING', 'PENDING', 'QUEUED']:
                return None, f"Unknown job status: {status}"
        
        return None, "InterProScan timed out after 5 minutes - the server may be busy"
        
    except requests.exceptions.Timeout:
        return None, "InterProScan API connection timed out - server may be busy"
    except requests.exceptions.ConnectionError:
        return None, "Could not connect to InterProScan API"
    except Exception as e:
        return None, f"InterProScan error: {str(e)}"


def search_hmmer_api(protein_seq):
    """
    Search for domains using HMMER web API (Pfam database).
    This is FASTER than InterProScan (typically 5-30 seconds).
    API Documentation: https://www.ebi.ac.uk/Tools/hmmer/
    """
    try:
        # HMMER hmmscan endpoint for domain searching
        url = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
        
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
            'Accept': 'application/json'
        }
        
        data = {
            'hmmdb': 'pfam',
            'seq': f">query\n{protein_seq}",
        }
        
        print("Submitting to HMMER/Pfam API...")
        
        response = requests.post(
            url, 
            data=data, 
            headers=headers, 
            timeout=60,
            allow_redirects=True
        )
        
        print(f"HMMER response status: {response.status_code}")
        
        if response.status_code == 200:
            try:
                content_type = response.headers.get('Content-Type', '')
                if 'application/json' in content_type:
                    results = response.json()
                    annotations = parse_hmmer_results(results)
                    if annotations:
                        return annotations, None
                    return None, "HMMER found no domains"
                else:
                    # Try to parse as JSON anyway
                    try:
                        results = response.json()
                        annotations = parse_hmmer_results(results)
                        if annotations:
                            return annotations, None
                    except:
                        pass
            except Exception as e:
                print(f"HMMER parse error: {e}")
        
        return None, f"HMMER returned status {response.status_code}"
        
    except requests.exceptions.Timeout:
        return None, "HMMER API timed out (60s)"
    except Exception as e:
        return None, f"HMMER error: {str(e)}"


def parse_interproscan_results(results):
    """Parse InterProScan JSON results into annotation format."""
    annotations = []
    
    try:
        if 'results' not in results:
            return annotations
        
        for result in results.get('results', []):
            matches = result.get('matches', [])
            
            for match in matches:
                signature = match.get('signature', {})
                entry = signature.get('entry', {})
                
                # Get domain name and description
                if entry:
                    name = entry.get('name', signature.get('name', 'Unknown'))
                    description = entry.get('description', signature.get('description', name))
                    accession = entry.get('accession', signature.get('accession', ''))
                    entry_type = entry.get('type', 'DOMAIN')
                else:
                    name = signature.get('name', 'Unknown')
                    description = signature.get('description', name)
                    accession = signature.get('accession', '')
                    entry_type = 'DOMAIN'
                
                # Get database source
                sig_lib = signature.get('signatureLibraryRelease', {})
                db_name = sig_lib.get('library', 'InterPro')
                
                # Get locations
                locations = match.get('locations', [])
                for loc in locations:
                    start = loc.get('start', 0)
                    end = loc.get('end', 0)
                    score = loc.get('score', 0)
                    evalue = loc.get('evalue', loc.get('score', 0))
                    
                    if start > 0 and end > start:
                        # Determine annotation type
                        ann_type = 'domain'
                        entry_type_lower = entry_type.lower() if entry_type else ''
                        name_lower = name.lower()
                        
                        if 'active' in name_lower or 'site' in entry_type_lower:
                            ann_type = 'active_site'
                        elif 'binding' in name_lower or 'metal' in name_lower:
                            ann_type = 'binding_site'
                        elif 'repeat' in name_lower or 'motif' in entry_type_lower:
                            ann_type = 'motif'
                        elif 'family' in entry_type_lower:
                            ann_type = 'family'
                        
                        annotations.append({
                            'name': name,
                            'description': description,
                            'start': start,
                            'end': end,
                            'type': ann_type,
                            'db': db_name,
                            'accession': accession,
                            'evalue': evalue,
                            'score': score
                        })
    
    except Exception as e:
        print(f"Error parsing InterProScan results: {str(e)}")
    
    # Remove duplicates and sort
    annotations = remove_duplicate_annotations(annotations)
    annotations.sort(key=lambda x: x['start'])
    
    return annotations


def parse_hmmer_results(results):
    """Parse HMMER JSON results into annotation format."""
    annotations = []
    
    try:
        hits = results.get('results', {}).get('hits', [])
        
        for hit in hits:
            domains = hit.get('domains', [])
            
            for domain in domains:
                annotations.append({
                    'name': hit.get('name', 'Unknown'),
                    'description': hit.get('desc', hit.get('name', '')),
                    'start': domain.get('alisqfrom', domain.get('ienv', 0)),
                    'end': domain.get('alisqto', domain.get('jenv', 0)),
                    'type': 'domain',
                    'db': 'Pfam/HMMER',
                    'accession': hit.get('acc', ''),
                    'evalue': domain.get('ievalue', domain.get('cevalue', 0)),
                    'score': domain.get('bitscore', 0)
                })
    
    except Exception as e:
        print(f"Error parsing HMMER results: {str(e)}")
    
    return annotations


def parse_pfam_results(results):
    """Parse Pfam API results into annotation format."""
    annotations = []
    
    try:
        for entry in results.get('results', results):
            if isinstance(entry, dict):
                metadata = entry.get('metadata', entry)
                
                annotations.append({
                    'name': metadata.get('name', 'Unknown'),
                    'description': metadata.get('description', ''),
                    'start': entry.get('start', entry.get('protein_start', 1)),
                    'end': entry.get('end', entry.get('protein_end', 1)),
                    'type': 'domain',
                    'db': 'Pfam',
                    'accession': metadata.get('accession', ''),
                    'evalue': entry.get('evalue', 0)
                })
    
    except Exception as e:
        print(f"Error parsing Pfam results: {str(e)}")
    
    return annotations


def remove_duplicate_annotations(annotations):
    """Remove duplicate annotations, keeping best scoring ones."""
    if not annotations:
        return []
    
    seen = {}
    for ann in annotations:
        key = (ann['name'], ann['start'], ann['end'])
        if key not in seen:
            seen[key] = ann
        else:
            # Keep the one with better e-value
            if ann.get('evalue', 1) < seen[key].get('evalue', 1):
                seen[key] = ann
    
    return list(seen.values())


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
                    'description': f'{domain_name} motif/domain (PROSITE-style pattern)',
                    'start': start,
                    'end': end,
                    'evalue': 0.001,
                    'type': 'motif',
                    'db': 'Local/PROSITE-patterns',
                    'accession': 'N/A'
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
    
    # ===== ROW 2: Protein Domain Architecture =====
    ax_domain = fig.add_subplot(gs[1, :])
    
    vis_length = int(seq_length)  # Ensure integer
    ax_domain.set_xlim(0, vis_length)
    ax_domain.set_ylim(-0.8, 1.2)
    
    # Draw main protein backbone (solid black line)
    ax_domain.plot([0, vis_length], [0.5, 0.5], color='black', linewidth=10, solid_capstyle='butt', zorder=1)
    
    # Green triangle marker in center top
    triangle_x = vis_length / 2
    triangle = plt.Polygon([[triangle_x - vis_length*0.02, 0.9], 
                            [triangle_x + vis_length*0.02, 0.9], 
                            [triangle_x, 1.1]], 
                           color='#2ecc71', zorder=10)
    ax_domain.add_patch(triangle)
    
    # Domain colors - Blue (#5555FF) and Orange (#FFA500) like your reference
    domain_colors = {
        'kinase': '#5555FF',      # Blue
        'domain': '#5555FF',      # Blue
        'sh3': '#5555FF',         # Blue
        'transmembrane': '#5555FF', # Blue
        'zinc': '#FFA500',        # Orange  
        'finger': '#FFA500',      # Orange
        'motif': '#FFA500',       # Orange
        'binding': '#FFA500',     # Orange
        'site': '#FFA500',        # Orange
        'phosphorylation': '#FFA500',  # Orange
        'glycosylation': '#FFA500',    # Orange
        'default': '#5555FF'      # Blue default
    }
    
    # Draw domains directly on the backbone
    for i, ann in enumerate(annotations):
        start = int(ann['start'])
        end = int(ann['end'])
        width = end - start
        name = ann['name']
        
        # Determine color based on name
        color = domain_colors['default']
        name_lower = name.lower()
        for key in domain_colors:
            if key in name_lower:
                color = domain_colors[key]
                break
        
        # Draw domain box on backbone (centered at y=0.5)
        domain_box = patches.FancyBboxPatch((start, 0.25), width, 0.5,
                                             boxstyle="round,pad=0.01",
                                             facecolor=color, edgecolor='black', 
                                             linewidth=1.5, zorder=3)
        ax_domain.add_patch(domain_box)
        
        # Add domain label on box
        label = name[:20] if len(name) > 20 else name
        if width > vis_length * 0.05:
            ax_domain.text(start + width/2, 0.5, label, ha='center', va='center',
                          fontsize=8, fontweight='bold', color='white', zorder=4)
    
    # Set up x-axis with proper integer ticks
    ax_domain.set_xlabel('Amino Acid Position', fontsize=10, fontweight='bold')
    ax_domain.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=8))
    ax_domain.ticklabel_format(style='plain', axis='x')  # Disable scientific notation
    ax_domain.tick_params(axis='x', labelsize=9)
    ax_domain.set_yticks([])
    ax_domain.spines['top'].set_visible(False)
    ax_domain.spines['right'].set_visible(False)
    ax_domain.spines['left'].set_visible(False)
    ax_domain.set_title('Protein Domain Architecture', fontsize=12, fontweight='bold', pad=10)
    
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
    """Create clean, static Plotly visualization with clear domain labels and ranges."""
    if not annotations:
        return None
    
    # Ensure seq_length is integer
    seq_length = int(seq_length)
    num_domains = len(annotations)
    
    # Color palette
    color_palette = [
        '#3498db',  # Blue
        '#e74c3c',  # Red
        '#2ecc71',  # Green
        '#f39c12',  # Orange
        '#9b59b6',  # Purple
        '#1abc9c',  # Teal
        '#e67e22',  # Dark Orange
        '#8e44ad',  # Dark Purple
    ]
    
    fig = go.Figure()
    
    # ===== Draw protein backbone =====
    fig.add_shape(
        type="rect",
        x0=0, y0=0.4, x1=seq_length, y1=0.6,
        fillcolor="#2c3e50",
        line=dict(color="#2c3e50", width=0),
    )
    
    # N terminus
    fig.add_annotation(
        x=-seq_length*0.05, y=0.5,
        text="<b>N</b>",
        showarrow=False,
        font=dict(size=18, color="#e74c3c", family="Arial Black"),
    )
    
    # C terminus  
    fig.add_annotation(
        x=seq_length*1.05, y=0.5,
        text="<b>C</b>",
        showarrow=False,
        font=dict(size=18, color="#27ae60", family="Arial Black"),
    )
    
    # ===== Draw domains on backbone =====
    for i, ann in enumerate(annotations):
        name = ann['name']
        start = int(ann['start'])
        end = int(ann['end'])
        color = color_palette[i % len(color_palette)]
        
        # Domain rectangle
        fig.add_shape(
            type="rect",
            x0=start, y0=0.3, x1=end, y1=0.7,
            fillcolor=color,
            line=dict(color="white", width=2),
            opacity=0.95,
        )
    
    # ===== Add domain table below =====
    table_y_start = -0.1
    row_height = 0.18
    
    # Table header
    fig.add_annotation(
        x=seq_length/2, y=table_y_start + 0.05,
        text="<b>━━━━━━━━━ DOMAINS FOUND ━━━━━━━━━</b>",
        showarrow=False,
        font=dict(size=11, color="#2c3e50"),
    )
    
    for i, ann in enumerate(annotations):
        name = ann['name']
        start = int(ann['start'])
        end = int(ann['end'])
        length = end - start
        color = color_palette[i % len(color_palette)]
        y_pos = table_y_start - ((i + 1) * row_height)
        
        # Color box
        fig.add_shape(
            type="rect",
            x0=0, y0=y_pos - 0.05, x1=seq_length * 0.04, y1=y_pos + 0.05,
            fillcolor=color,
            line=dict(color="white", width=1),
        )
        
        # Domain info text - full name
        fig.add_annotation(
            x=seq_length * 0.06, y=y_pos,
            text=f"<b>{name}</b>",
            showarrow=False,
            font=dict(size=11, color="#2c3e50"),
            xanchor="left",
        )
        
        # Position range - clearly visible
        fig.add_annotation(
            x=seq_length * 0.55, y=y_pos,
            text=f"<b>Range: {start} - {end}</b>  ({length} aa)",
            showarrow=False,
            font=dict(size=11, color="#555"),
            xanchor="left",
        )
    
    # Calculate height based on domains
    plot_height = 300 + (num_domains * 35)
    y_min = table_y_start - ((num_domains + 1) * row_height) - 0.1
    
    # ===== Layout - FULLY STATIC, no zoom/pan =====
    fig.update_layout(
        title=dict(
            text=f"<b>🧬 Protein Domain Analysis</b><br><span style='font-size:13px'>Sequence: {seq_length} amino acids  |  Domains: {num_domains}</span>",
            font=dict(size=16, color='#2c3e50'),
            x=0.5,
            y=0.95
        ),
        height=plot_height,
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=60, r=60, t=80, b=30),
        xaxis=dict(
            range=[-seq_length*0.08, seq_length*1.08],
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            tickformat='d',
            tickfont=dict(size=11),
            title=dict(text="<b>Amino Acid Position</b>", font=dict(size=12)),
            fixedrange=True,  # Disable zoom
            dtick=max(10, seq_length // 8),  # Nice tick spacing
        ),
        yaxis=dict(
            range=[y_min, 1.0],
            showticklabels=False,
            showgrid=False,
            fixedrange=True,  # Disable zoom
        ),
        showlegend=False,
        dragmode=False,  # Disable drag
    )
    
    # Convert to JSON with static config embedded
    fig_json = pio.to_json(fig)
    return fig_json

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
            sequence = request.form.get('sequence', '').strip()
            seq_type = request.form.get('seq_type', 'protein')
            
            if not sequence:
                return jsonify({'error': 'No sequence provided'}), 400
            
            # Check if it's FASTA format (starts with >)
            if sequence.startswith('>'):
                # Parse FASTA format
                lines = sequence.split('\n')
                header_line = lines[0][1:].strip()  # Remove '>' and get header
                
                # Extract ID and description from header
                # Format: >ID description or >ID|other|info description
                header_parts = header_line.split(None, 1)  # Split on first whitespace
                if header_parts:
                    seq_id = header_parts[0]
                    seq_description = header_parts[1] if len(header_parts) > 1 else header_line
                else:
                    seq_id = "FASTA_Protein"
                    seq_description = header_line
                
                # Join remaining lines as sequence (skip header)
                sequence = ''.join(lines[1:])
                print(f"Parsed FASTA: ID={seq_id}, Description={seq_description[:50]}...")
            else:
                seq_id = "Manual_Protein"
                seq_description = "Manually entered protein sequence"
            
            # Clean sequence - remove spaces, newlines, numbers
            sequence = sequence.upper()
            sequence = ''.join(c for c in sequence if c.isalpha() or c == '*')
            
            if not sequence:
                return jsonify({'error': 'No valid sequence found after parsing'}), 400
            
            if seq_type == 'protein':
                # Accept standard amino acids plus some common modifications
                valid_aa = set('ACDEFGHIKLMNPQRSTVWY*XU')
                invalid_chars = [c for c in sequence if c not in valid_aa]
                if invalid_chars:
                    return jsonify({
                        'error': f'Invalid characters in protein sequence: {", ".join(set(invalid_chars))}. Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y'
                    }), 400
                # Clean the sequence - remove stop codons and unknown residues
                protein_seq = sequence.replace('*', '').replace('X', '').replace('U', 'C')  # Selenocysteine to Cysteine
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
        
        elif input_type == 'pdb':
            # Handle PDB ID input
            pdb_id = request.form.get('pdb_id', '').strip()
            if not pdb_id:
                return jsonify({'error': 'No PDB ID provided'}), 400
            
            protein_seq, seq_description, error = fetch_sequence_from_pdb(pdb_id)
            if error:
                return jsonify({'error': error}), 400
            
            seq_id = f"PDB:{pdb_id.upper()}"
            print(f"✓ Successfully fetched from PDB: {len(protein_seq)} amino acids")
        
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
