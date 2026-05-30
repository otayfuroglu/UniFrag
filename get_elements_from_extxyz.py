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
        
    sorted_elements = sorted(list(unique_elements))
    
    # Save elements to text file
    print(f"\nFound {len(sorted_elements)} unique element(s): {', '.join(sorted_elements)}")
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            for element in sorted_elements:
                f.write(f"{element}\n")
        print(f"Successfully wrote elements list to '{output_path}'.")
    except Exception as e:
        print(f"Error writing to output file '{output_path}': {e}")


if __name__ == "__main__":
    main()
