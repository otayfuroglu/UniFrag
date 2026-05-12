# UniFrag

UniFrag is a Python framework designed to cleanly and robustly extract chemically-meaningful cluster models from Metal-Organic Frameworks (MOFs), Covalent Organic Frameworks (COFs), and Biological Macromolecules (PDBs). It is specifically optimized for creating stable models for Density Functional Theory (DFT) and QM calculations.

MOF2frag is a Python script designed to cleanly and robustly extract chemically-meaningful symmetric cluster models from Metal-Organic Framework (MOF) `.cif` files. It is specifically optimized for creating stable models for Density Functional Theory (DFT) calculations by using an intelligently-anchored Graph Breadth-First-Search.

## Features
- **Symmetric Cluster Extraction**: Targets a core metal cluster (e.g., an Mg3 or Ni3 node) and strictly completes all organic ligands bonded to those core metals. Disconnected or grazing organic linkers are gracefully ignored to prevent floating geometries.
- **Intact Organic Linkers**: The bond-graph traversal guarantees that your ligands are pulled completely symmetrically without slicing through covalent bonds.
- **Automatic Valence Capping**: Any carboxylate, phenolic, or bound Oxygen atom that loses at least one metal during extraction is automatically "capped" with an artificial Hydrogen atom placed precisely along the averaged bond vector at 0.96 Å. This creates standard and perfectly neutralized closed-shell clusters automatically!

## Installation
Dependencies:
```bash
pip install pymatgen numpy openbabel
conda install -c conda-forge pdbfixer openmm
```

## Usage
Simply run the Python script on your CIF file.
```bash
python fragmenter.py <path_to_cif> --radius 4.0 --output my_fragment.xyz
```

### Arguments
- `cif_path`: Path to your MOF `.cif` file.
- `--center`: The atomic index in the unit cell that acts as the physical origin. Leave it blank (or set to `-1`) and the script will automatically lock onto the most centralized metal atom to avoid periodic boundary wrapping.
- `--radius`: The spatial sphere (in Ångströms) defining the central "core metals". **Radius = 4.0 Å is highly recommended** for MOF-74 analogues to consistently isolate precisely 3 core metal atoms. *Warning: Radii larger than 4.5 Å will trigger a structural warning as they may disrupt molecular symmetries depending on your specific MOF lattice.*
- `--output`: Filepath to save the resulting `.xyz` fragment.

## Example
If you want to extract a 4-linker, 3-Metal core out of an MOF-74 variant:
```bash
python fragmenter.py example/MgMOF74_clean_fromCORE.cif --radius 4.0 --output Mg_cluster.xyz
```
This produces an `.xyz` trajectory containing 3 Mg atoms, 32 C atoms, 24 O atoms, and 23 capping/native Hydrogens!

## Biological Macromolecules (BioMolFragmenter)

UniFrag now supports the fragmentation of single-chain biological macromolecules (proteins) using a sequence-based sliding window approach. This is ideal for extracting contiguous peptide segments for QM calculations.

**Key Features:**
- **Sliding Window:** Iterates along the backbone sequence using `--window-size` and `--stride` to produce well-defined fragments without relying on spatial radii.
- **Strict Neutralization:** Biological fragments are strictly neutralized for QM calculations. The internal bond graph automatically detects charged functional groups at pH 7 (e.g., carboxylates, ammonium, guanidinium, imidazolium) and adds or removes single Hydrogen atoms to enforce a formal charge of exactly **0**.
- **Cut-Bond Only Capping:** Heavy atoms are taken directly from the PDB coordinates and are never modified. Only the severed peptide bonds at the N- and C-terminal cuts are capped with single Hydrogen atoms.
- **PDBFixer Integration:** Automatically runs `pdbfixer` to resolve missing residues, missing heavy atoms, and missing native Hydrogens at a physiological pH before fragmentation begins.
- **Chemical Deduplication:** Uses a strict chemical fingerprint to detect and skip writing redundant `.xyz` files when highly repetitive sequences produce identical geometries.

### Usage
Run the fragmentation script specifying `--kind bio`. You can process individual PDB files or entire directories using the shell runner.

```bash
# Process a single PDB file with window size 5, stride 1, and pH 7.0
python fragmentation_oop.py example_protein.pdb --kind bio --window-size 5 --stride 1 --ph 7.0

# Process an entire directory of PDBs using the shell helper
./run_bio_family.sh test_on_bio_mol 5 1 7.0
```

### Bio Arguments
- `--window-size`: Number of amino acid residues to include in each fragment (default: 5). **Decrease this number to reduce the total number of atoms in your fragments.**
- `--stride`: The step size in residues between consecutive windows (default: 1).
- `--ph`: The physiological pH to use when PDBFixer adds missing hydrogens (default: 7.0).
- `--no-pdbfixer`: Skip the PDBFixer pre-processing step and use the original raw PDB file directly.
