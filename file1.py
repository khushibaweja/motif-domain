from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
import json
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import plotly.graph_objects as go
import plotly.offline as pyo
import os
import sys

# Set your email for NCBI compliance (Replace with your actual email)
Entrez.email = os.environ.get('NCBI_EMAIL', 'your.email@example.com')
if Entrez.email == 'your.email@example.com':
    print("⚠️  Warning: Using default email. Set NCBI_EMAIL environment variable for better compliance.")

# Dynamic separator width
SEPARATOR_WIDTH = 80

def display_welcome_message():
    """Display a welcome message with program options."""
    print("\n" + "=" * SEPARATOR_WIDTH)
    print("🧬 PROTEIN DOMAIN VISUALIZATION TOOL")
    print("=" * SEPARATOR_WIDTH)
    print("This tool helps you analyze protein sequences and visualize domains/motifs.")
    print("\nChoose your input method:")

def get_user_choice():
    """Get user input for sequence source with detailed options."""
    print("\n" + "═" * SEPARATOR_WIDTH)
    print("SELECT INPUT METHOD")
    print("═" * SEPARATOR_WIDTH)
    print("1. 📁 Upload FASTA file")
    print("2. 🔍 Fetch by NCBI Accession Number") 
    print("3. ✍  Enter sequence manually")
    print("4. 🚪 Exit program")
    
    while True:
        try:
            choice = input("\nEnter your choice (1-4): ").strip()
            if choice in ['1', '2', '3', '4']:
                return choice
            else:
                print("❌ Invalid choice. Please enter 1, 2, 3, or 4.")
        except KeyboardInterrupt:
            print("\n\n👋 Program interrupted by user. Goodbye!")
            sys.exit(0)

def upload_fasta_file():
    """Handle FASTA file upload with error checking."""
    print("\n" + "─" * SEPARATOR_WIDTH)
    print("📁 FASTA FILE UPLOAD")
    print("─" * SEPARATOR_WIDTH)
    print("Note: Both protein and nucleotide FASTA files are supported.")
    
    while True:
        file_path = input("Enter the path to your FASTA file: ").strip().strip('"')
        
        if not file_path:
            print("❌ No file path provided. Please try again.")
            continue
            
        if not os.path.exists(file_path):
            print(f"❌ File not found: {file_path}")
            retry = input("Would you like to try another path? (y/n): ").lower()
            if retry != 'y':
                return None
            continue
            
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            if not records:
                print("❌ No valid sequences found in the file.")
                retry = input("Would you like to try another file? (y/n): ").lower()
                if retry != 'y':
                    return None
                continue
                
            print(f"✅ Successfully loaded {len(records)} sequence(s) from {os.path.basename(file_path)}")
            
            # Check sequence type
            for record in records:
                seq_str = str(record.seq).upper()
                if all(aa in 'ACDEFGHIKLMNPQRSTVWY*' for aa in seq_str):
                    print(f"   - {record.id}: Protein sequence ({len(record.seq)} aa)")
                else:
                    print(f"   - {record.id}: Nucleotide sequence ({len(record.seq)} bp)")
            
            return records
            
        except Exception as e:
            print(f"❌ Error reading FASTA file: {e}")
            retry = input("Would you like to try another file? (y/n): ").lower()
            if retry != 'y':
                return None

def fetch_by_accession():
    """Fetch sequence by NCBI accession number."""
    print("\n" + "─" * SEPARATOR_WIDTH)
    print("🔍 FETCH BY ACCESSION NUMBER")
    print("─" * SEPARATOR_WIDTH)
    print("Supported: Protein (NP_123456, P12345) and Nucleotide (NM_123456) accessions")
    
    while True:
        accession = input("\nEnter accession number: ").strip()
        
        if not accession:
            print("❌ No accession number provided.")
            retry = input("Try again? (y/n): ").lower()
            if retry != 'y':
                return None
            continue
            
        print(f"🔍 Searching for {accession}...")
        record = fetch_sequence_from_accession(accession)
        
        if record:
            return [record]
        else:
            retry = input("❌ Failed to fetch sequence. Try another accession? (y/n): ").lower()
            if retry != 'y':
                return None

def enter_sequence_manually():
    """Allow user to enter sequence manually."""
    print("\n" + "─" * SEPARATOR_WIDTH)
    print("✍  MANUAL SEQUENCE ENTRY")
    print("─" * SEPARATOR_WIDTH)
    print("Choose sequence type:")
    print("1. Protein sequence (one-letter amino acid codes)")
    print("2. DNA sequence (nucleotide bases)")
    
    while True:
        seq_type = input("\nEnter choice (1-2): ").strip()
        if seq_type in ['1', '2']:
            break
        else:
            print("❌ Invalid choice. Please enter 1 or 2.")
    
    if seq_type == '1':
        print("\nEnter protein sequence using one-letter codes (ACDEFGHIKLMNPQRSTVWY):")
        print("Example: MGHHHHHHSSGVDLGTENLYFQSQ")
    else:
        print("\nEnter DNA sequence (ACGTN bases):")
        print("Example: ATGGGCCAGATTGTC")
    
    sequences = []
    
    while True:
        print("\nEnter sequence (or 'done' to finish, 'cancel' to exit):")
        sequence_input = input().strip().upper()
        
        if sequence_input == 'DONE':
            if sequences:
                break
            else:
                print("❌ No sequences entered. Please enter at least one sequence.")
                continue
        elif sequence_input == 'CANCEL':
            return None
        elif not sequence_input:
            print("❌ Empty sequence. Please enter a valid sequence.")
            continue
            
        # Validate sequence based on type
        if seq_type == '1':
            # Protein sequence validation
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY*')
            if all(aa in valid_aa for aa in sequence_input):
                # Create a mock record
                record = SeqRecord(
                    Seq(sequence_input.replace('*', '')),
                    id=f"manual_protein_{len(sequences)+1}",
                    description="Manually entered protein sequence"
                )
                sequences.append(record)
                print(f"✅ Protein sequence {len(sequences)} added ({len(sequence_input)} amino acids)")
                
                another = input("Add another sequence? (y/n): ").lower()
                if another != 'y':
                    break
            else:
                invalid_chars = set(sequence_input) - valid_aa
                print(f"❌ Invalid amino acid characters: {', '.join(invalid_chars)}")
                print("Valid characters: ACDEFGHIKLMNPQRSTVWY*")
        else:
            # DNA sequence validation
            valid_bases = set('ACGTN')
            if all(base in valid_bases for base in sequence_input):
                record = SeqRecord(
                    Seq(sequence_input),
                    id=f"manual_dna_{len(sequences)+1}",
                    description="Manually entered DNA sequence"
                )
                sequences.append(record)
                print(f"✅ DNA sequence {len(sequences)} added ({len(sequence_input)} bases)")
                
                another = input("Add another sequence? (y/n): ").lower()
                if another != 'y':
                    break
            else:
                invalid_chars = set(sequence_input) - valid_bases
                print(f"❌ Invalid nucleotide characters: {', '.join(invalid_chars)}")
                print("Valid characters: ACGTN")
    
    return sequences

def safe_translate(dna_seq):
    """Translate DNA sequence to protein safely."""
    try:
        # Remove any non-coding characters and ensure length is multiple of 3
        clean_seq = ''.join([base for base in str(dna_seq).upper() if base in 'ACGTN'])
        
        if len(clean_seq) < 3:
            print("⚠️  DNA sequence too short for translation (minimum 3 bases required)")
            return None
            
        remainder = len(clean_seq) % 3
        if remainder != 0:
            print(f"⚠️  Trimming {remainder} bases to make sequence length a multiple of 3")
            clean_seq = clean_seq[:-remainder]  # Trim to multiple of 3
        
        protein_seq = Seq(clean_seq).translate(to_stop=True)
        
        if len(protein_seq) == 0:
            print("⚠️  Translation resulted in empty protein (stop codon at start)")
            return None
            
        return protein_seq
    except Exception as e:
        print(f"❌ Error in translation: {e}")
        return None

def analyze_protein_features(protein_seq):
    """Analyze basic protein features."""
    print("\n" + "=" * SEPARATOR_WIDTH)
    print("PROTEIN FEATURE ANALYSIS")
    print("=" * SEPARATOR_WIDTH)
    seq = str(protein_seq).upper()
    
    if len(seq) == 0:
        print("Protein sequence is empty.")
        return
    
    # Amino acid counts
    aa_counts = {}
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        aa_counts[aa] = seq.count(aa)
    
    # Feature calculations
    charged = aa_counts['R'] + aa_counts['K'] + aa_counts['D'] + aa_counts['E']
    hydrophobic = aa_counts['A'] + aa_counts['V'] + aa_counts['L'] + aa_counts['I'] + aa_counts['M'] + aa_counts['F'] + aa_counts['W']
    polar = aa_counts['S'] + aa_counts['T'] + aa_counts['N'] + aa_counts['Q']
    special = aa_counts['C'] + aa_counts['G'] + aa_counts['P'] + aa_counts['H'] + aa_counts['Y']
    
    length = len(seq)
    
    print(f"Sequence length: {length} amino acids")
    print(f"Charged residues (R,K,D,E): {charged} ({charged/length*100:.1f}%)")
    print(f"Hydrophobic residues: {hydrophobic} ({hydrophobic/length*100:.1f}%)")
    print(f"Polar residues: {polar} ({polar/length*100:.1f}%)")
    print(f"Special residues: {special} ({special/length*100:.1f}%)")
    
    # Molecular weight approximation (average amino acid MW ~110 Da)
    approx_mw = length * 110
    print(f"Approximate molecular weight: {approx_mw:.0f} Da")
    
    # Basic predictions
    if hydrophobic/length > 0.4:
        print("Note: High hydrophobicity - possible transmembrane regions")
    if aa_counts['C'] >= 2:
        print(f"Note: Contains {aa_counts['C']} cysteine residues - potential disulfide bonds")

def search_motifs_interproscan(protein_seq):
    """Search for motifs and domains using InterProScan API or mock data."""
    print("\n" + "=" * SEPARATOR_WIDTH)
    print("Domain and Motif Search")
    print("=" * SEPARATOR_WIDTH)
    
    seq_length = len(protein_seq)
    
    if seq_length < 20:
        print("Sequence is too short for domain analysis.")
        return []
    
    # Use mock data for demonstration
    print("\n⚠️  WARNING: Using mock domain data for demonstration purposes")
    print("Real protein domain analysis requires InterProScan API setup or other bioinformatics tools.")
    print("These are simulated domains for visualization demonstration only.\n")
    
    # Generate mock domains based on sequence length
    mock_annotations = []
    
    # Add some common protein domains as mock data
    possible_domains = [
        {'name': 'Kinase domain', 'type': 'domain', 'fraction': 0.3},
        {'name': 'SH3 domain', 'type': 'domain', 'fraction': 0.1},
        {'name': 'Zinc finger', 'type': 'motif', 'fraction': 0.05},
        {'name': 'Transmembrane region', 'type': 'domain', 'fraction': 0.2},
        {'name': 'Active site', 'type': 'active_site', 'fraction': 0.02},
    ]
    
    current_pos = 1
    for domain in possible_domains:
        if current_pos >= seq_length * 0.8:  # Don't go beyond 80% of sequence
            break
            
        domain_length = int(seq_length * domain['fraction'])
        if current_pos + domain_length > seq_length:
            domain_length = seq_length - current_pos - 10
            
        if domain_length > 10:  # Only add domains of reasonable size
            end_pos = current_pos + domain_length
            mock_annotations.append({
                'name': domain['name'],
                'description': f"Mock {domain['type']}",
                'start': current_pos,
                'end': end_pos,
                'type': domain['type'],
                'db': 'MockDB',
                'evalue': '1e-10',
                'score': '100'
            })
            current_pos = end_pos + int(seq_length * 0.05)  # Add gap between domains
    
    if mock_annotations:
        print(f"Found {len(mock_annotations)} domains/motifs:")
        for i, ann in enumerate(mock_annotations, 1):
            print(f"  {i}. {ann['name']} ({ann['type']})")
            print(f"     Positions: {ann['start']}-{ann['end']}")
    else:
        print("No domains or motifs found.")
    
    return mock_annotations

def create_ncbi_style_visualization(seq_length, annotations, title="Conserved Domain Visualization"):
    """Create NCBI CDD-style domain visualization."""
    
    # Ensure seq_length is reasonable (proteins are rarely > 5000 aa)
    seq_length = int(seq_length)
    if seq_length > 50000:
        print(f"Warning: Very long sequence ({seq_length}). Might be nucleotide length?")
        # If it looks like nucleotide length, divide by 3
        if seq_length > 100000:
            seq_length = seq_length // 3
            print(f"Adjusted to protein length: {seq_length} aa")
    
    vis_length = max(seq_length, 100)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8), 
                                   gridspec_kw={'height_ratios': [1, 1.5]})
    
    fig.patch.set_facecolor('white')

    # --- Top panel: Domain architecture ---
    ax1.set_xlim(0, vis_length)
    ax1.set_ylim(0, 3)
    
    # Draw sequence line (black backbone)
    ax1.plot([0, vis_length], [2.0, 2.0], color='black', linewidth=6)
    
    # Add green triangle marker at center
    center_x = vis_length / 2
    triangle = plt.Polygon([[center_x - vis_length*0.02, 2.5], 
                           [center_x + vis_length*0.02, 2.5], 
                           [center_x, 2.8]], 
                          color='#2ecc71', zorder=10)
    ax1.add_patch(triangle)

    # Draw domains and motifs - Blue (#5555FF) and Orange (#FFA500)
    for ann in annotations:
        start = int(ann['start'])
        end = int(ann['end'])
        width = end - start
        ann_name = ann['name'].lower()
        
        # Color coding - blue for domains, orange for motifs/binding sites
        if any(keyword in ann_name for keyword in ['motif', 'site', 'zinc', 'finger', 'binding', 'phosph', 'glyco']):
            color = '#FFA500'  # Orange
        else:
            color = '#5555FF'  # Blue
        
        # Draw domain box on backbone
        rect = patches.FancyBboxPatch((start, 1.6), width, 0.8,
                                     boxstyle="round,pad=0.01",
                                     facecolor=color, edgecolor='black', 
                                     linewidth=1.5, zorder=3)
        ax1.add_patch(rect)
        
        # Add label for larger domains
        if width > vis_length * 0.05:
            ax1.text(start + width/2, 2.0, ann['name'][:15], 
                    ha='center', va='center', fontsize=8, fontweight='bold', color='white')

    # Add position markers with proper integer formatting
    ax1.set_xlabel('Amino Acid Position', fontweight='bold')
    ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=8))
    ax1.ticklabel_format(style='plain', axis='x')  # Disable scientific notation
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_title("Protein Domain Architecture", fontsize=14, fontweight='bold', pad=20)

    # --- Bottom panel: Simple features track ---
    ax2.set_xlim(0, vis_length)
    ax2.set_ylim(0, len(annotations) + 1)
    
    for i, ann in enumerate(annotations):
        y_pos = len(annotations) - i
        start = int(ann['start'])
        end = int(ann['end'])
        ann_name = ann['name'].lower()
        
        # Color based on type - blue for domains, orange for motifs
        if any(keyword in ann_name for keyword in ['motif', 'site', 'zinc', 'finger', 'binding', 'phosph', 'glyco']):
            color = '#FFA500'  # Orange
        else:
            color = '#5555FF'  # Blue
        
        rect = patches.FancyBboxPatch((start, y_pos - 0.4), end-start, 0.8,
                                     boxstyle="round,pad=0.01",
                                     facecolor=color, alpha=0.9, edgecolor='black')
        ax2.add_patch(rect)
        ax2.text(start + (end-start)/2, y_pos, ann['name'][:20], 
                ha='center', va='center', fontsize=8, color='white', fontweight='bold')

    ax2.set_xlabel('Amino Acid Position', fontweight='bold')
    ax2.set_ylabel('Features', fontweight='bold')
    ax2.set_title('Domain and Motif Features', fontsize=12, fontweight='bold')
    ax2.set_yticks(range(1, len(annotations) + 1))
    ax2.set_yticklabels([ann['name'][:25] for ann in annotations])
    ax2.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=8))
    ax2.ticklabel_format(style='plain', axis='x')  # Disable scientific notation
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    return fig

def create_interactive_track_plot(seq_length, annotations, title="Interactive Motif and Domain Tracks"):
    """
    Creates an INTERACTIVE track-style plot using Plotly.
    Each annotation is on its own track with hover-over details.
    """
    if not annotations:
        print("No annotations to plot interactively.")
        return None

    fig = go.Figure()
    
    # Color map
    color_map = {
        'domain': 'blue',
        'motif': 'red',
        'active_site': 'red',
        'metal_binding': 'green',
        'other': 'gray'
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

        # Determine type for coloring
        ann_type = ann.get('type', 'domain')
        if any(keyword in name.lower() for keyword in ['active', 'site', 'catalytic']):
            ann_type = 'active_site'
        elif any(keyword in name.lower() for keyword in ['metal', 'binding', 'zinc']):
            ann_type = 'metal_binding'
        elif any(keyword in name.lower() for keyword in ['phosphorylation', 'motif']):
            ann_type = 'motif'
        elif 'domain' in name.lower():
            ann_type = 'domain'
        
        color = color_map.get(ann_type, 'other')
        
        # Create the hover text
        hover_text = (
            f"<b>{ann['name']}</b><br>"
            f"Type: {ann_type}<br>"
            f"Position: {start} - {end}<br>"
            f"Length: {end - start} aa<br>"
            f"Source: {ann.get('db', 'N/A')}<br>"
            f"Description: {ann.get('description', 'N/A')}"
        )

        # Add the domain as a shape
        fig.add_shape(
            type="rect",
            x0=start, y0=y_pos - 0.4,
            x1=end, y1=y_pos + 0.4,
            line=dict(color="Black", width=1),
            fillcolor=color,
            opacity=0.8
        )
        
        # Add transparent bar for hover events
        fig.add_trace(go.Bar(
            x=[(end + start) / 2],
            y=[y_pos],
            width=(end - start),
            base=start,
            orientation='h',
            marker_color=color,
            opacity=0.0,
            text=hover_text,
            hoverinfo='text',
            name=name
        ))

    # Configure layout
    fig.update_layout(
        title=f'<b>{title}</b>',
        xaxis_title="Amino Acid Position",
        yaxis_title="Domains and Motifs",
        height=max(600, len(annotations) * 40),
        xaxis=dict(range=[0, seq_length]),
        yaxis=dict(
            tickmode='array',
            tickvals=y_positions,
            ticktext=y_labels
        ),
        showlegend=False,
        plot_bgcolor='white',
        xaxis_gridcolor='lightgray',
        xaxis_gridwidth=1,
    )
    
    # Save to HTML file
    plot_filename = f"{title.replace(' ', '_')}_interactive.html"
    pyo.plot(fig, filename=plot_filename, auto_open=True)
    print(f"Interactive plot saved as: {plot_filename}")
    
    return fig

def visualize_annotations(seq_length, annotations, title="Domain Visualization"):
    """Main visualization function."""
    if not annotations:
        print("⚠️  No annotations to visualize.")
        return
        
    print("\nGenerating visualizations...")
    
    try:
        # Ask which visualizations to create
        print("\nVisualization options:")
        print("1. Static plot only")
        print("2. Interactive plot only")
        print("3. Both (default)")
        
        viz_choice = input("Choose visualization type (1-3, or Enter for both): ").strip()
        
        if viz_choice == '1' or viz_choice == '3' or viz_choice == '':
            # Create static visualization
            print("Generating static visualization...")
            fig1 = create_ncbi_style_visualization(seq_length, annotations, title)
            filename = f"{title.replace(' ', '_')}.png"
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"✅ Static plot saved as: {filename}")
            plt.show()
        
        if viz_choice == '2' or viz_choice == '3' or viz_choice == '':
            # Create interactive visualization
            print("Generating interactive visualization...")
            create_interactive_track_plot(seq_length, annotations, title)
    except Exception as e:
        print(f"❌ Error creating visualizations: {e}")

def fetch_sequence_from_accession(acc):
    """Fetch a sequence from NCBI using accession number."""
    try:
        print(f"🔍 Fetching sequence {acc} from NCBI...")
        
        # First try protein database
        try:
            handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            print(f"✅ Successfully fetched PROTEIN: {record.description}")
            return record
        except:
            # If protein fails, try nucleotide database
            print("🔄 Trying nucleotide database...")
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            print(f"✅ Successfully fetched NUCLEOTIDE: {record.description}")
            return record
            
    except Exception as e:
        print(f"❌ Error fetching accession {acc}: {e}")
        return None

def process_sequence_record(record, sequence_number):
    """Process a single sequence record."""
    print(f"\n{'=' * (SEPARATOR_WIDTH + 20)}")
    print(f"ANALYZING SEQUENCE {sequence_number}: {record.id}")
    print(f"{'=' * (SEPARATOR_WIDTH + 20)}")
    print(f"Description: {record.description}")
    
    # Check if it's DNA or protein
    seq = record.seq
    seq_str = str(seq).upper()
    is_protein = all(aa in 'ACDEFGHIKLMNPQRSTVWY*' for aa in seq_str)
    
    if is_protein:
        print(f"Type: Protein sequence")
        print(f"Length: {len(seq)} amino acids")
        protein_str = str(seq).replace('*', '')
        analyze_protein_features(protein_str)
        
        # Use mock data for motif and domain finding
        annotations = search_motifs_interproscan(protein_str)
    else:
        print(f"Type: DNA sequence")
        print(f"Length: {len(seq)} nucleotides")
        
        # Translate DNA to protein
        protein_seq = safe_translate(seq)
        if protein_seq:
            protein_str = str(protein_seq)
            print(f"Translated protein length: {len(protein_str)} amino acids")
            analyze_protein_features(protein_str)
            
            # Use mock data for motif and domain finding
            annotations = search_motifs_interproscan(protein_str)
        else:
            print("❌ Failed to translate DNA sequence to protein.")
            return
    
    print("\n" + "=" * SEPARATOR_WIDTH)
    print("GENERATING VISUALIZATIONS")
    print("=" * SEPARATOR_WIDTH)
    
    if annotations:
        create_viz = input("Generate visualizations? (y/n): ").lower()
        if create_viz == 'y':
            visualize_annotations(len(protein_str), annotations, 
                                 f"Domains in {record.id}")
        else:
            print("Skipping visualizations.")
    else:
        print("No annotations to visualize.")

def main():
    """Main function with enhanced user interface."""
    try:
        display_welcome_message()
    except Exception as e:
        print(f"❌ Error displaying welcome message: {e}")
        return
    
    while True:
        try:
            choice = get_user_choice()
        except Exception as e:
            print(f"❌ Error getting user choice: {e}")
            continue
        
        query_records = None
        try:
            if choice == '1':
                query_records = upload_fasta_file()
            elif choice == '2':
                query_records = fetch_by_accession()
            elif choice == '3':
                query_records = enter_sequence_manually()
            elif choice == '4':
                print("\n👋 Thank you for using the Protein Domain Visualization Tool!")
                sys.exit(0)
        except KeyboardInterrupt:
            print("\n\n⚠️  Operation cancelled by user.")
            continue
        except Exception as e:
            print(f"\n❌ Error processing input: {e}")
            continue
        
        if not query_records:
            print("\n❌ No sequences to process.")
            continue
            
        # Process all sequences
        for i, record in enumerate(query_records, 1):
            try:
                process_sequence_record(record, i)
                
                # Ask if user wants to continue with next sequence
                if i < len(query_records):
                    cont = input(f"\nProcess next sequence? ({i+1}/{len(query_records)}) (y/n): ").lower()
                    if cont != 'y':
                        break
                        
            except Exception as e:
                print(f"❌ Error processing sequence {i}: {e}")
                continue
        
        # Ask if user wants to analyze more sequences
        another = input("\n🔁 Would you like to analyze more sequences? (y/n): ").lower()
        if another != 'y':
            print("\n👋 Thank you for using the Protein Domain Visualization Tool!")
            break

if __name__ == "__main__":
    main()