import os
import glob
import shutil
import warnings
import networkx as nx
from pymatgen.core import Structure
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN

# Set environment variable for CCDC
os.environ["CSD_DATA_DIRECTORY"] = "/Users/omert/CCDC/ccdc-data/csd"
from ccdc.io import EntryReader

# Ignore warnings from pymatgen
warnings.filterwarnings("ignore")

def clean_guests():
    cif_dir = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/zn_cr_cifs_noduplicated/cifs"
    backup_dir = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/zn_cr_cifs_noduplicated/cifs_backup_guests"
    report_path = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/guest_removal_report.md"
    brain_report_path = "/Users/omert/.gemini/antigravity/brain/153d2da2-7e4a-4474-b36a-6be8db573d0d/guest_removal_report.md"
    
    if not os.path.exists(backup_dir):
        os.makedirs(backup_dir)
        
    cif_paths = sorted(glob.glob(os.path.join(cif_dir, "*.cif")))
    print(f"Starting guest detection and removal on {len(cif_paths)} purified parent CIFs...")
    
    cleaned_records = []
    failed_records = []
    
    for path in cif_paths:
        refcode = os.path.splitext(os.path.basename(path))[0]
        backup_path = os.path.join(backup_dir, f"{refcode}.cif")
        
        # If backup exists, it means this structure has already been cleaned!
        # We can construct the record from the backup (original) and current (cleaned) files.
        if os.path.exists(backup_path):
            try:
                # Read original (backup)
                reader_orig = EntryReader(backup_path)
                entry_orig = reader_orig[0]
                components_orig = entry_orig.molecule.components
                guest_comps = []
                for i, comp in enumerate(components_orig):
                    has_zn = any(a.atomic_symbol == 'Zn' for a in comp.atoms)
                    if not has_zn:
                        guest_comps.append(comp.formula)
                
                struct_orig = Structure.from_file(backup_path, occupancy_tolerance=100.0)
                orig_formula = struct_orig.formula
                
                # Read cleaned (current)
                struct_clean = Structure.from_file(path, occupancy_tolerance=100.0)
                clean_formula = struct_clean.formula
                
                # Verification with CCDC
                reader_verify = EntryReader(path)
                entry_verify = reader_verify[0]
                verify_components = entry_verify.molecule.components
                remaining_guests = []
                for i, comp in enumerate(verify_components):
                    has_zn = any(a.atomic_symbol == 'Zn' for a in comp.atoms)
                    if not has_zn:
                        remaining_guests.append(comp.formula)
                        
                num_removed_sites = len(struct_orig) - len(struct_clean)
                
                # Check status
                if remaining_guests:
                    status = "Warning (Guests remain)"
                else:
                    status = "Success"
                    
                cleaned_records.append({
                    "refcode": refcode,
                    "original_formula": orig_formula,
                    "cleaned_formula": clean_formula,
                    "removed_guests": ", ".join(guest_comps),
                    "num_removed_sites": num_removed_sites,
                    "status": status
                })
                continue
            except Exception as e:
                print(f"  Error loading backup/current structure for {refcode}: {str(e)}")
                failed_records.append({
                    "refcode": refcode,
                    "error": f"Failed to load backup/cleaned file: {str(e)}"
                })
                continue
                
        # Normal flow if backup doesn't exist
        try:
            # 1. CCDC Pre-filtering
            reader = EntryReader(path)
            entry = reader[0]
            mol = entry.molecule
            components = mol.components
            
            guest_comps = []
            for i, comp in enumerate(components):
                has_zn = any(a.atomic_symbol == 'Zn' for a in comp.atoms)
                if not has_zn:
                    guest_comps.append(comp.formula)
            
            if not guest_comps:
                # No guests, skip
                continue
                
            print(f"\n>>> REFCODE {refcode} has guest molecules: {', '.join(guest_comps)}")
            
            # 2. Backup the original file
            shutil.copy(path, backup_path)
            print(f"  Backed up original to {backup_path}")
            
            # 3. Read with Pymatgen
            struct = Structure.from_file(path, occupancy_tolerance=100.0)
            orig_formula = struct.formula
            
            # Verify Zn presence
            if not any(site.species_string == 'Zn' for site in struct):
                raise ValueError("No Zinc atoms found in the structure")
                
            indices_to_remove = set()
            method_used = None
            
            try:
                # Build periodic structure graph with JmolNN
                from pymatgen.analysis.local_env import JmolNN
                sg_jmol = StructureGraph.with_local_env_strategy(struct, JmolNN())
                g_jmol = sg_jmol.graph.to_undirected()
                pmg_components_jmol = list(nx.connected_components(g_jmol))
                
                indices_to_keep = []
                for comp in pmg_components_jmol:
                    has_zn = any(struct[idx].species_string == 'Zn' for idx in comp)
                    if has_zn:
                        indices_to_keep.extend(comp)
                        
                indices_to_remove = set(range(len(struct))) - set(indices_to_keep)
                method_used = "JmolNN"
                
            except Exception as e_jmol:
                print(f"  JmolNN failed: {str(e_jmol)}. Falling back to CrystalNN...")
                indices_to_remove = set()
                method_used = None
                
            # Fallback to CrystalNN
            if method_used is None:
                try:
                    cnn = CrystalNN(search_cutoff=3.0)
                    sg = StructureGraph.with_local_env_strategy(struct, cnn)
                    g = sg.graph.to_undirected()
                    pmg_components = list(nx.connected_components(g))
                    
                    indices_to_keep = []
                    for comp in pmg_components:
                        has_zn = any(struct[idx].species_string == 'Zn' for idx in comp)
                        if has_zn:
                            indices_to_keep.extend(comp)
                            
                    indices_to_remove = set(range(len(struct))) - set(indices_to_keep)
                    method_used = "CrystalNN"
                    print(f"  Successfully fell back to CrystalNN for component isolation.")
                except Exception as e_cnn:
                    print(f"  CrystalNN fallback also failed: {str(e_cnn)}")
                    indices_to_remove = set()
            
            if not indices_to_remove:
                print(f"  Warning: No guest components isolated in either CrystalNN or JmolNN for {refcode}.")
                status = "Warning (Not cleaned)"
                cleaned_records.append({
                    "refcode": refcode,
                    "original_formula": orig_formula,
                    "cleaned_formula": orig_formula,
                    "removed_guests": ", ".join(guest_comps),
                    "num_removed_sites": 0,
                    "status": status
                })
                continue
                
            print(f"  Removing {len(indices_to_remove)} sites from {len(struct)} total sites using {method_used}...")
            
            # Remove sites
            clean_struct = struct.copy()
            clean_struct.remove_sites(sorted(list(indices_to_remove), reverse=True))
            
            # Overwrite original CIF in-place
            clean_struct.to(filename=path)
            print(f"  Cleaned CIF written successfully. New formula: {clean_struct.formula}")
            
            # 4. Post-removal verification with CCDC
            reader_verify = EntryReader(path)
            entry_verify = reader_verify[0]
            verify_components = entry_verify.molecule.components
            
            remaining_guests = []
            for i, comp in enumerate(verify_components):
                has_zn = any(a.atomic_symbol == 'Zn' for a in comp.atoms)
                if not has_zn:
                    remaining_guests.append(comp.formula)
                    
            if remaining_guests:
                print(f"  Verification WARNING: {refcode} still has guest components: {', '.join(remaining_guests)}")
                status = f"Warning (Guests remain, using {method_used})"
            else:
                print(f"  Verification SUCCESS: {refcode} is fully purified.")
                status = f"Success (using {method_used})"
                
            cleaned_records.append({
                "refcode": refcode,
                "original_formula": orig_formula,
                "cleaned_formula": clean_struct.formula,
                "removed_guests": ", ".join(guest_comps),
                "num_removed_sites": len(indices_to_remove),
                "status": status
            })
            
        except Exception as e:
            print(f"  ERROR processing {refcode}: {str(e)}")
            failed_records.append({
                "refcode": refcode,
                "error": str(e)
            })
            
    # Write summary reports
    print(f"\nPurification complete. Cleaned {len(cleaned_records)} structures. Failed on {len(failed_records)} structures.")
    
    # Generate Report Content
    report_content = f"# Guest Molecule Removal Report\n\n"
    report_content += f"**Total parent CIFs scanned:** {len(cif_paths)}  \n"
    report_content += f"**Purified Zn MOFs requiring cleaning:** {len(cleaned_records) + len(failed_records)}  \n"
    report_content += f"**Successfully cleaned and verified:** {len([r for r in cleaned_records if 'Success' in r['status']])}  \n"
    report_content += f"**Cleaned with warnings/failures:** {len([r for r in cleaned_records if 'Success' not in r['status']]) + len(failed_records)}  \n\n"
    
    if cleaned_records:
        # Sort by status and refcode
        cleaned_records = sorted(cleaned_records, key=lambda x: (x["status"], x["refcode"]))
        report_content += "## Cleaned parent CIFs\n\n"
        report_content += "| REFCODE | Original Formula | Cleaned Formula | Removed Guests | Removed Atoms | Status |\n"
        report_content += "|---|---|---|---|---|---|\n"
        for r in cleaned_records:
            report_content += f"| {r['refcode']} | {r['original_formula']} | {r['cleaned_formula']} | {r['removed_guests']} | {r['num_removed_sites']} | {r['status']} |\n"
        report_content += "\n"
        
    if failed_records:
        report_content += "## Failed/Error cases\n\n"
        report_content += "| REFCODE | Error Message |\n"
        report_content += "|---|---|\n"
        for r in failed_records:
            report_content += f"| {r['refcode']} | {r['error']} |\n"
        report_content += "\n"
        
    # Save reports
    with open(report_path, "w") as f:
        f.write(report_content)
    with open(brain_report_path, "w") as f:
        f.write(report_content)
    print(f"Reports saved to:\n  - {report_path}\n  - {brain_report_path}")

if __name__ == "__main__":
    clean_guests()
