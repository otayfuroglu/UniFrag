#!/usr/bin/env python3
"""
Utility script to extract unique element types from an ExtXYZ file and write them to a text file.

Requirements:
    - ase (Atomic Simulation Environment)
"""

import argparse
from pathlib import Path

try:
    from ase.io import read
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

# Periodic Table atomic numbers mapping for element sorting
ATOMIC_NUMBERS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
    'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}


def get_atomic_number(symbol):
    """Retrieves the atomic number of an element symbol, falling back to 999 if unknown."""
    normalized = symbol.strip().capitalize()
    return ATOMIC_NUMBERS.get(normalized, 999)


def main():
    parser = argparse.ArgumentParser(
        description="Extract unique chemical element types from an ExtXYZ file."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        default="fragments_collection.extxyz",
        help="Path to the input ExtXYZ file (default: fragments_collection.extxyz)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="elements.txt",
        help="Path to the output text file (default: elements.txt)"
    )
    
    args = parser.parse_args()
    
    if not ASE_AVAILABLE:
        print("Error: The 'ase' package is required to run this script.")
        print("Please install it using: pip install ase")
        return
        
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file '{input_path}' does not exist.")
        return
        
    output_path = Path(args.output)
    
    print(f"Reading structures from '{input_path}'...")
    try:
        # Read all frames in the file
        frames = read(str(input_path), index=":")
    except Exception as e:
        print(f"Error reading ExtXYZ file: {e}")
        return
        
    n_frames = len(frames)
    print(f"Loaded {n_frames} structure(s).")
    if n_frames == 0:
        print("No structures found in input file.")
        return
        
    # Extract unique elements
    unique_elements = set()
    for frame in frames:
        unique_elements.update(frame.get_chemical_symbols())
        
    # Sort elements by atomic number (smallest atomic number first)
    sorted_elements = sorted(list(unique_elements), key=get_atomic_number)
    
    # Save elements to text file
    print(f"\nFound {len(sorted_elements)} unique element(s) (sorted by atomic number): {', '.join(sorted_elements)}")
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            for element in sorted_elements:
                f.write(f"{element}\n")
        print(f"Successfully wrote elements list to '{output_path}'.")
    except Exception as e:
        print(f"Error writing to output file '{output_path}': {e}")


if __name__ == "__main__":
    main()

