#!/usr/bin/env python3
"""
Utility script to separate structures in an ExtXYZ file into two categories:
1. Structures with <= N atoms (smaller or equal)
2. Structures with > N atoms (larger)

Where N is a threshold value specified by the user.

Requirements:
    - ase (Atomic Simulation Environment)
"""

import argparse
from pathlib import Path

try:
    from ase.io import read, write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


def main():
    parser = argparse.ArgumentParser(
        description="Separate structures inside an ExtXYZ file into exactly two files based on an atom count threshold."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        default="fragments_collection.extxyz",
        help="Path to the input ExtXYZ file (default: fragments_collection.extxyz)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default="separated_by_atoms",
        help="Directory to save the separated ExtXYZ files (default: separated_by_atoms)"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=int,
        required=True,
        help="Atom count threshold value N. Structures with <= N atoms are put in the smaller file; structures with > N atoms in the larger file."
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
        
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
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
        
    threshold = args.threshold
    smaller_frames = []
    larger_frames = []
    
    # Categorize structures based on the threshold
    for frame in frames:
        n_atoms = len(frame)
        if n_atoms <= threshold:
            smaller_frames.append(frame)
        else:
            larger_frames.append(frame)
            
    # File paths
    smaller_file = output_dir / f"smaller_or_equal_to_{threshold}.extxyz"
    larger_file = output_dir / f"larger_than_{threshold}.extxyz"
    
    # Save the files
    print(f"\nCategorizing structures based on threshold N = {threshold}:")
    
    if smaller_frames:
        try:
            write(str(smaller_file), smaller_frames)
            print(f"  - {smaller_file.name:<32} : Saved {len(smaller_frames):>4} structure(s) (<= {threshold} atoms)")
        except Exception as e:
            print(f"  - Error saving {smaller_file.name}: {e}")
    else:
        print(f"  - No structures found with <= {threshold} atoms (no file written).")
        
    if larger_frames:
        try:
            write(str(larger_file), larger_frames)
            print(f"  - {larger_file.name:<32} : Saved {len(larger_frames):>4} structure(s) (> {threshold} atoms)")
        except Exception as e:
            print(f"  - Error saving {larger_file.name}: {e}")
    else:
        print(f"  - No structures found with > {threshold} atoms (no file written).")
        
    print("\nPost-processing complete!")


if __name__ == "__main__":
    main()
