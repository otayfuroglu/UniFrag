import os
import sys
import glob
import re
import shutil
import argparse
import warnings
import concurrent.futures
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN

# Ignore warnings from pymatgen
warnings.filterwarnings("ignore")

# Standard metals list from the project (as fallback / reference)
METALS_SET = {
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
    'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Al', 'Ga', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po'
}

def is_metal_element(symbol):
    """
    Checks if an element symbol represents a metal.
    Uses Pymatgen Element properties with a robust fallback.
    """
    try:
        el = Element(symbol)
        if el.is_metal or el.is_transition_metal or el.is_lanthanoid or el.is_actinoid:
            return True
    except Exception:
        pass
    return symbol in METALS_SET

def get_base_refcode(filename):
    """
    Extracts the 6-character base CSD refcode from a filename.
    """
    name = os.path.splitext(os.path.basename(filename))[0].upper()
    tokens = re.split(r'[_ \-]+', name)
    for token in tokens:
        match = re.search(r'([A-Z]{6})', token)
        if match:
            return match.group(1)
    
    alpha_only = re.sub(r'[^A-Z]', '', name)
    if len(alpha_only) >= 6:
        return alpha_only[:6]
    return name

def filter_single_cif(args_tuple):
    """
    Phase 1 Worker: Parses CIF, filters target metal, mixed metals, and excluded elements.
    """
    path, metal, exclude_elements, single_metal_only = args_tuple
    filename = os.path.basename(path)
    refcode = get_base_refcode(filename)
    try:
        struct = Structure.from_file(path, occupancy_tolerance=100.0)
    except Exception as e:
        return {"status": "error", "filename": filename, "error": str(e)}
        
    all_elements = {el.symbol for el in struct.composition.elements}
    if metal not in all_elements:
        return {"status": "no_target_metal", "filename": filename}
        
    if single_metal_only:
        metals_present = {el for el in all_elements if is_metal_element(el)}
        if len(metals_present) > 1:
            return {"status": "mixed_metals", "filename": filename, "metals": sorted(list(metals_present))}
            
    excluded_present = all_elements.intersection(set(exclude_elements))
    if excluded_present:
        return {"status": "excluded_elements", "filename": filename, "excluded": sorted(list(excluded_present))}
        
    return {
        "status": "valid",
        "path": path,
        "filename": filename,
        "refcode": refcode
    }

def make_ordered_structure(struct):
    """
    Returns an ordered version of the structure by taking the highest-occupancy
    specie for any disordered/multi-species sites.
    """
    if struct.is_ordered:
        return struct
        
    from pymatgen.core import PeriodicSite
    new_sites = []
    for site in struct:
        if not site.is_ordered:
            best_specie = max(site.species.items(), key=lambda x: x[1])[0]
            new_sites.append(PeriodicSite(best_specie, site.frac_coords, site.lattice, properties=site.properties))
        else:
            new_sites.append(site)
            
    return Structure.from_sites(new_sites)

def process_guest_stripping(args_tuple):
    """
    Phase 3 Worker: Performs connected component analysis to remove guests/solvent molecules.
    """
    path, filename, refcode, dest_dir, metal, remove_guests, deduplicate_refcodes = args_tuple
    dest_filename = f"{refcode}.cif" if deduplicate_refcodes else filename
    dest_path = os.path.join(dest_dir, dest_filename)
    backup_dir = os.path.join(dest_dir, "cifs_backup_guests")
    
    try:
        struct = Structure.from_file(path, occupancy_tolerance=100.0)
        orig_formula = struct.formula
    except Exception as e:
        return {
            "status": "error",
            "filename": filename,
            "error": f"Failed to re-parse structure: {str(e)}"
        }
        
    if not remove_guests:
        shutil.copy(path, dest_path)
        return {"status": "skipped", "filename": filename}
        
    # Create ordered copy for graph analysis to avoid Pymatgen specie errors on disordered sites
    ordered_struct = make_ordered_structure(struct)
    
    indices_to_remove = set()
    method_used = None
    e_cnn = None
    
    # Try CrystalNN first
    try:
        import networkx as nx
        cnn = CrystalNN(search_cutoff=3.0)
        sg = StructureGraph.with_local_env_strategy(ordered_struct, cnn)
        g = sg.graph.to_undirected()
        pmg_components = list(nx.connected_components(g))
        
        indices_to_keep = []
        for comp in pmg_components:
            has_metal = any(ordered_struct[idx].species_string == metal for idx in comp)
            if has_metal:
                indices_to_keep.extend(comp)
                
        indices_to_remove = set(range(len(struct))) - set(indices_to_keep)
        if indices_to_remove:
            method_used = "CrystalNN"
    except Exception as ex:
        e_cnn = ex
        indices_to_remove = set()
        
    # Fallback to JmolNN
    if not indices_to_remove:
        try:
            import networkx as nx
            from pymatgen.analysis.local_env import JmolNN
            sg_jmol = StructureGraph.with_local_env_strategy(ordered_struct, JmolNN())
            g_jmol = sg_jmol.graph.to_undirected()
            pmg_components_jmol = list(nx.connected_components(g_jmol))
            
            indices_to_keep = []
            for comp in pmg_components_jmol:
                has_metal = any(ordered_struct[idx].species_string == metal for idx in comp)
                if has_metal:
                    indices_to_keep.extend(comp)
                    
            indices_to_remove = set(range(len(struct))) - set(indices_to_keep)
            if indices_to_remove:
                method_used = "JmolNN"
        except Exception as e_jmol:
            err_msg = f"CrystalNN failed ({str(e_cnn)}) and JmolNN failed ({str(e_jmol)})"
            shutil.copy(path, dest_path)
            return {
                "status": "warning",
                "filename": filename,
                "refcode": refcode,
                "error": err_msg
            }
            
    if indices_to_remove and method_used:
        # Guests found and removed!
        os.makedirs(backup_dir, exist_ok=True)
        backup_path = os.path.join(backup_dir, filename)
        shutil.copy(path, backup_path)
        
        clean_struct = struct.copy()
        clean_struct.remove_sites(sorted(list(indices_to_remove), reverse=True))
        clean_struct.to(filename=dest_path)
        
        return {
            "status": "cleaned",
            "filename": filename,
            "refcode": refcode,
            "original_formula": orig_formula,
            "cleaned_formula": clean_struct.formula,
            "removed_atoms": len(indices_to_remove),
            "method": method_used
        }
    else:
        # No guests found
        shutil.copy(path, dest_path)
        return {"status": "no_guests", "filename": filename}

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    parser = argparse.ArgumentParser(description="Standard Pre-processing Pipeline before Fragmentation")
    parser.add_argument("--src_dir", required=True, help="Source directory containing raw CIF files")
    parser.add_argument("--dest_dir", required=True, help="Destination directory for clean pre-processed CIFs")
    parser.add_argument("--metal", default="Zn", help="Target metal element symbol (default: Zn)")
    parser.add_argument("--exclude_elements", nargs="+", default=['I', 'Si', 'Br', 'Se', 'As'],
                        help="List of element symbols to filter out/exclude")
    parser.add_argument("--single_metal_only", type=str2bool, nargs='?', const=True, default=True,
                        help="Skip mixed-metal MOFs containing other metals (default: True)")
    parser.add_argument("--remove_guests", type=str2bool, nargs='?', const=True, default=True,
                        help="Detect and strip guest/solvent/counter-ion molecules (default: True)")
    parser.add_argument("--deduplicate_refcodes", type=str2bool, nargs='?', const=True, default=True,
                        help="Group by 6-letter parent refcode and keep only the best representative (default: True)")
    parser.add_argument("--brain_dir", default=None, help="Optional brain artifacts directory to copy the report to")
    
    args = parser.parse_args()
    
    os.makedirs(args.dest_dir, exist_ok=True)
    cif_paths = sorted(glob.glob(os.path.join(args.src_dir, "*.cif")))
    
    if not cif_paths:
        print(f"Error: No CIF files found in source directory {args.src_dir}")
        sys.exit(1)
        
    print(f"Pre-processing Pipeline initialized:")
    print(f"  Source dir: {args.src_dir}")
    print(f"  Dest dir: {args.dest_dir}")
    print(f"  Target metal: {args.metal}")
    print(f"  Excluded elements: {args.exclude_elements}")
    print(f"  Single metal only: {args.single_metal_only}")
    print(f"  Remove guests: {args.remove_guests}")
    print(f"  Deduplicate refcodes: {args.deduplicate_refcodes}")
    print(f"  Total files to scan: {len(cif_paths)}")
    print("-" * 60)
    
    # ----------------------------------------------------
    # Phase 1: Parallel Filtering
    # ----------------------------------------------------
    filter_args = [(p, args.metal, args.exclude_elements, args.single_metal_only) for p in cif_paths]
    
    parsing_errors = {}
    no_target_metal = []
    mixed_metals = {}
    excluded_elements_found = {}
    valid_structures = []
    
    print("\nPhase 1: Filtering structures (parallel execution)...")
    num_cpus = os.cpu_count() or 2
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpus) as executor:
        results = list(executor.map(filter_single_cif, filter_args))
        
    for res in results:
        status = res["status"]
        fn = res["filename"]
        if status == "error":
            parsing_errors[fn] = res["error"]
        elif status == "no_target_metal":
            no_target_metal.append(fn)
        elif status == "mixed_metals":
            mixed_metals[fn] = res["metals"]
        elif status == "excluded_elements":
            excluded_elements_found[fn] = res["excluded"]
        elif status == "valid":
            valid_structures.append(res)
            
    print(f"Phase 1 Complete:")
    print(f"  Parsing errors: {len(parsing_errors)}")
    print(f"  Missing target metal ({args.metal}): {len(no_target_metal)}")
    print(f"  Mixed metal skipped: {len(mixed_metals)}")
    print(f"  Excluded elements skipped: {len(excluded_elements_found)}")
    print(f"  Valid structures remaining: {len(valid_structures)}")
    
    # ----------------------------------------------------
    # Phase 2: Deduplication
    # ----------------------------------------------------
    selected_structures = []
    deduplicated_groups = {}
    
    if args.deduplicate_refcodes:
        print("\nPhase 2: Deduplicating base refcodes...")
        groups = {}
        for item in valid_structures:
            groups.setdefault(item["refcode"], []).append(item)
            
        for refcode, items in sorted(groups.items()):
            # Priority: ASR (0) > FSR (1) > Ion (2) > others (3), then alphabetical
            def sort_key(x):
                fn = x["filename"].upper()
                if '_ASR' in fn:
                    pri = 0
                elif '_FSR' in fn:
                    pri = 1
                elif '_ION' in fn:
                    pri = 2
                else:
                    pri = 3
                return (pri, x["filename"])
                
            sorted_items = sorted(items, key=sort_key)
            selected = sorted_items[0]
            selected_structures.append(selected)
            
            skipped = [x["filename"] for x in sorted_items[1:]]
            if skipped:
                deduplicated_groups[refcode] = {
                    "selected": selected["filename"],
                    "skipped": skipped
                }
        print(f"Phase 2 Complete:")
        print(f"  Unique parent refcodes: {len(selected_structures)}")
        print(f"  Duplicate files skipped: {sum(len(v['skipped']) for v in deduplicated_groups.values())}")
    else:
        selected_structures = valid_structures
        print("\nPhase 2: Deduplication skipped per user choice.")
        
    # ----------------------------------------------------
    # Phase 3: Guest Stripping & Connected Component Analysis
    # ----------------------------------------------------
    print(f"\nPhase 3: Connected Component analysis & Guest Stripping (parallel execution on {len(selected_structures)} files)...")
    guest_args = [
        (item["path"], item["filename"], item["refcode"], args.dest_dir, args.metal, args.remove_guests, args.deduplicate_refcodes)
        for item in selected_structures
    ]
    
    guest_cleaned_records = []
    guest_failed_records = []
    no_guests_count = 0
    skipped_count = 0
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpus) as executor:
        guest_results = list(executor.map(process_guest_stripping, guest_args))
        
    for res in guest_results:
        status = res["status"]
        if status == "cleaned":
            guest_cleaned_records.append(res)
        elif status == "error" or status == "warning":
            guest_failed_records.append(res)
        elif status == "no_guests":
            no_guests_count += 1
        elif status == "skipped":
            skipped_count += 1
            
    print(f"Phase 3 Complete:")
    print(f"  Successfully cleaned structures: {len(guest_cleaned_records)}")
    print(f"  Structures with no guests: {no_guests_count}")
    print(f"  Failed/Warnings: {len(guest_failed_records)}")
    
    # ----------------------------------------------------
    # Report Generation
    # ----------------------------------------------------
    print("\nGenerating Pre-processing Report...")
    report_path = os.path.join(args.dest_dir, "preprocess_report.md")
    
    report_content = f"# Pre-processing Pipeline Summary Report\n\n"
    report_content += f"This report was automatically generated by the UniFrag pre-processing pipeline script (`preprocess_dataset.py`).\n\n"
    
    report_content += "## Pipeline Parameters\n\n"
    report_content += f"- **Source Directory:** `{args.src_dir}`\n"
    report_content += f"- **Destination Directory:** `{args.dest_dir}`\n"
    report_content += f"- **Target Metal:** `{args.metal}`\n"
    report_content += f"- **Excluded Elements:** `{', '.join(args.exclude_elements)}`\n"
    report_content += f"- **Single Metal Only:** `{args.single_metal_only}`\n"
    report_content += f"- **Remove Guests:** `{args.remove_guests}`\n"
    report_content += f"- **Deduplicate Refcodes:** `{args.deduplicate_refcodes}`\n\n"
    
    report_content += "## Execution Statistics Summary\n\n"
    report_content += f"- **Total Crystal Structures Scanned:** {len(cif_paths)}\n"
    report_content += f"- **Failed to Parse:** {len(parsing_errors)}\n"
    report_content += f"- **No Target Metal ({args.metal}):** {len(no_target_metal)}\n"
    report_content += f"- **Mixed-Metal Excluded:** {len(mixed_metals)}\n"
    report_content += f"- **Excluded Elements Excluded:** {len(excluded_elements_found)}\n"
    if args.deduplicate_refcodes:
        report_content += f"- **Conformer/Duplicate Files Skipped:** {sum(len(v['skipped']) for v in deduplicated_groups.values())}\n"
    report_content += f"- **Final Output Structures Saved:** {len(selected_structures)}\n"
    if args.remove_guests:
        report_content += f"  - **Cleaned (Guests Detected & Stripped):** {len(guest_cleaned_records)}\n"
        report_content += f"  - **No Guests Detected (Intact):** {no_guests_count}\n"
        report_content += f"  - **Guest removal failed/warning:** {len(guest_failed_records)}\n"
    report_content += "\n"
    
    # 1. Parsing Errors Table
    if parsing_errors:
        report_content += "## Parsing Errors\n\n"
        report_content += "| Filename | Error Message |\n"
        report_content += "| :--- | :--- |\n"
        for fn, err in sorted(parsing_errors.items()):
            report_content += f"| {fn} | {err} |\n"
        report_content += "\n"
        
    # 2. Missing Target Metal List
    if no_target_metal:
        report_content += f"## Missing Target Metal ({args.metal})\n\n"
        report_content += f"Total structures skipped: **{len(no_target_metal)}**\n\n"
        report_content += "<details>\n<summary>Click to view files list</summary>\n\n"
        for fn in sorted(no_target_metal):
            report_content += f"- {fn}\n"
        report_content += "\n</details>\n\n"
        
    # 3. Mixed Metals Table
    if mixed_metals:
        report_content += "## Mixed Metals Excluded\n\n"
        report_content += "| Filename | Metals Found |\n"
        report_content += "| :--- | :--- |\n"
        for fn, metals in sorted(mixed_metals.items()):
            report_content += f"| {fn} | {', '.join(metals)} |\n"
        report_content += "\n"
        
    # 4. Excluded Elements Table
    if excluded_elements_found:
        report_content += "## Excluded Elements Found & Skipped\n\n"
        report_content += "| Filename | Unwanted Elements |\n"
        report_content += "| :--- | :--- |\n"
        for fn, elems in sorted(excluded_elements_found.items()):
            report_content += f"| {fn} | {', '.join(elems)} |\n"
        report_content += "\n"
        
    # 5. Deduplication Detail Table
    if deduplicated_groups:
        report_content += "## Deduplication & Conformer Filtering Detail\n\n"
        report_content += "| Parent REFCODE | Selected File | Skipped Duplicate Files |\n"
        report_content += "| :--- | :--- | :--- |\n"
        for refcode, details in sorted(deduplicated_groups.items()):
            report_content += f"| **{refcode}** | {details['selected']} | {', '.join(details['skipped'])} |\n"
        report_content += "\n"
        
    # 6. Guest Removal Table
    if guest_cleaned_records:
        report_content += "## Guest Molecule Removal Detail\n\n"
        report_content += "| REFCODE | Original File | Original Formula | Cleaned Formula | Removed Atoms | Method |\n"
        report_content += "| :--- | :--- | :--- | :--- | :--- | :--- |\n"
        for r in sorted(guest_cleaned_records, key=lambda x: x["refcode"]):
            report_content += f"| **{r['refcode']}** | {r['filename']} | {r['original_formula']} | {r['cleaned_formula']} | {r['removed_atoms']} | {r['method']} |\n"
        report_content += "\n"
        
    # 7. Guest Failures Table
    if guest_failed_records:
        report_content += "## Guest Stripping Warnings / Errors\n\n"
        report_content += "| Filename | Warning / Error Details |\n"
        report_content += "| :--- | :--- |\n"
        for r in sorted(guest_failed_records, key=lambda x: x["filename"]):
            report_content += f"| {r['filename']} | {r.get('error', 'Unknown error')} |\n"
        report_content += "\n"
        
    with open(report_path, "w") as f:
        f.write(report_content)
        
    # Also save a copy in the brain artifacts directory if specified or if fallback exists
    brain_dir = args.brain_dir
    if not brain_dir:
        # Fallback to the default path for backwards compatibility
        brain_dir = "/Users/omert/.gemini/antigravity/brain/153d2da2-7e4a-4474-b36a-6be8db573d0d"
    if brain_dir and os.path.exists(brain_dir):
        try:
            shutil.copy(report_path, os.path.join(brain_dir, "preprocess_report.md"))
        except Exception:
            pass
        
    print(f"\nReport written to {report_path}")
    print("-" * 60)
    print("Pre-processing pipeline execution completed successfully!")

if __name__ == "__main__":
    main()
