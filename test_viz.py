"""Test the visualization dashboard functions"""
import sys
sys.path.insert(0, '.')

# Import the functions
from app import calculate_protein_properties, create_visualization_static

# Test sequence
test_seq = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIA'
test_annotations = [
    {'name': 'Kinase Domain', 'start': 10, 'end': 40, 'type': 'Domain', 'description': 'Protein kinase domain'},
    {'name': 'ATP binding', 'start': 25, 'end': 35, 'type': 'Motif', 'description': 'ATP binding site'}
]

# Test properties calculation
print("Testing calculate_protein_properties...")
props = calculate_protein_properties(test_seq)
print(f"  MW: {props['mw']:.1f} kDa")
print(f"  pI: {props['pi']:.2f}")
print(f"  Cys count: {props['cys_count']}")
print("  OK!")

# Test visualization
print("\nTesting create_visualization_static (dashboard)...")
result = create_visualization_static(len(test_seq), test_annotations, test_seq, 'TEST_ID')
print(f"  Dashboard generated: {len(result)} characters")
print("  OK!")

print("\n*** SUCCESS - All visualization functions working! ***")
