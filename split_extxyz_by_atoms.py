#!/usr/bin/env python3
"""
Utility script to separate structures in an ExtXYZ file by their atom counts.

This script parses a multi-frame ExtXYZ file (such as fragments_collection.extxyz)
and groups the structures either by their exact number of atoms or by configurable
atom count ranges. The grouped structures are saved into separate ExtXYZ files.

Requirements:
    - ase (Atomic Simulation Environment)
"""

import os
import argparse
from pathlib import Path
from collections import defaultdict

try:
    from ase.io import read, write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


def parse_ranges(range_str):
    """Parses a comma-separated string of integers into a sorted list of bounds."""
    if not range_str:
        return []
    try:
        bounds = sorted(list(set(int(x.strip()) for x in range_str.split(",") if x.strip())))
        if any(b <= 0 for b in bounds):
            raise ValueError("All bounds must be positive integers.")
        return bounds
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Invalid range bounds list: {e}")


def group_by_exact(frames):
    """Groups frames by their exact number of atoms."""
    grouped = defaultdict(list)
    for frame in frames:
        n_atoms = len(frame)
        grouped[n_atoms].append(frame)
    return grouped


def group_by_range(frames, bounds):
    """Groups frames by atom count ranges based on bounds.
    
    e.g., bounds = [50, 100, 200]
    Ranges will be:
      - 0 to 50
      - 51 to 100
      - 101 to 200
      - 201 and above
    """
    grouped = defaultdict(list)
    
    for frame in frames:
        n_atoms = len(frame)
        placed = False
        
        # Find which range it belongs to
        for i, bound in enumerate(bounds):
            lower = bounds[i-1] + 1 if i > 0 else 0
            if lower <= n_atoms <= bound:
                range_key = f"{lower}_to_{bound}"
                grouped[range_key].append(frame)
                placed = True
                break
        
        if not placed:
            lower = bounds[-1] + 1 if bounds else 0
            range_key = f"{lower}_and_above"
            grouped[range_key].append(frame)
            
    return grouped


def main():
    parser = argparse.ArgumentParser(
        description="Separate structures inside an ExtXYZ file by atom count."
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
        "-m", "--mode",
        choices=["exact", "range"],
        default="exact",
        help="Separation mode: 'exact' count files or 'range' binned files (default: exact)"
    )
    parser.add_argument(
        "-r", "--ranges",
        type=parse_ranges,
        default="",
        help="Comma-separated range upper bounds (e.g., '50,100,200') for range mode."
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
        
    # Group frames based on chosen mode
    if args.mode == "exact":
        print("Grouping structures by exact atom count...")
        grouped = group_by_exact(frames)
    else:
        # If range mode is selected but no ranges were provided, default to a sensible set
        bounds = args.ranges
        if not bounds:
            bounds = [50, 100, 200, 500]
            print(f"No ranges provided; defaulting to bounds: {bounds}")
        print(f"Grouping structures by atom count ranges with bounds: {bounds}...")
        grouped = group_by_range(frames, bounds)
        
    # Write output files
    print(f"\nWriting separated structures to '{output_dir}/':")
    for key, frame_list in sorted(grouped.items(), key=lambda x: (isinstance(x[0], int), x[0])):
        count = len(frame_list)
        if args.mode == "exact":
            filename = f"atoms_{key}.extxyz"
            desc = f"Exactly {key} atoms"
        else:
            filename = f"range_{key}.extxyz"
            desc = f"Range {key} atoms"
            
        file_path = output_dir / filename
        try:
            write(str(file_path), frame_list)
            print(f"  - {filename:<25} : Saved {count:>4} structure(s) ({desc})")
        except Exception as e:
            print(f"  - Error saving {filename}: {e}")
            
    print("\nPost-processing complete!")


if __name__ == "__main__":
    main()
