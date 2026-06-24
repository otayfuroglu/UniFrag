import os
import glob
import pandas as pd
os.environ["CSD_DATA_DIRECTORY"] = "/Users/omert/CCDC/ccdc-data/csd"
from ccdc.io import EntryReader

def main():
    cif_dir = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/zn_cr_cifs_noduplicated/cifs"
    output_csv_path = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/zn_low_coordination_parents.csv"
    artifact_csv_path = "/Users/omert/.gemini/antigravity/brain/153d2da2-7e4a-4474-b36a-6be8db573d0d/zn_low_coordination_parents.csv"
    
    cif_paths = sorted(glob.glob(os.path.join(cif_dir, "*.cif")))
    print(f"Scanning {len(cif_paths)} parent CIFs for Zn coordination numbers 0, 1, and 3...")
    
    rows = []
    
    for idx, path in enumerate(cif_paths):
        refcode = os.path.splitext(os.path.basename(path))[0].upper()
        
        try:
            reader = EntryReader(path)
            entry = reader[0]
            ccdc_mol = entry.molecule
            
            low_coord_envs = []
            for atom in ccdc_mol.atoms:
                if atom.atomic_symbol == 'Zn':
                    neighbors = [n.atomic_symbol for n in atom.neighbours if n.atomic_symbol != 'H']
                    cn = len(neighbors)
                    if cn in [0, 1, 3]:
                        counts = {el: neighbors.count(el) for el in set(neighbors)}
                        formula = "".join(f"{el}{counts[el]}" for el in sorted(counts.keys()))
                        key = f"CN={cn} [{formula}]" if cn > 0 else "CN=0 [None]"
                        low_coord_envs.append(key)
                        
            if low_coord_envs:
                # Deduplicate environment strings for the report
                env_str = "; ".join(sorted(list(set(low_coord_envs))))
                rows.append({
                    'cif_file': os.path.basename(path),
                    'refcode': refcode,
                    'low_coordination_environments': env_str
                })
        except Exception as e:
            print(f"Error processing {refcode}: {e}")
            
    # Write to CSV
    df = pd.DataFrame(rows)
    df.to_csv(output_csv_path, index=False)
    df.to_csv(artifact_csv_path, index=False)
    
    print(f"\nDone! Identified {len(df)} parent MOFs containing Zn with CN in [0, 1, 3].")
    print(f"Results saved to:")
    print(f"  - {output_csv_path}")
    print(f"  - {artifact_csv_path}")

if __name__ == "__main__":
    main()
