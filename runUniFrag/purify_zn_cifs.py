import os
import glob
import shutil
os.environ["CSD_DATA_DIRECTORY"] = "/Users/omert/CCDC/ccdc-data/csd"
from ccdc.io import EntryReader

def main():
    cif_dir = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/zn_cr_cifs_noduplicated/cifs"
    target_dir = "/Users/omert/Desktop/UniFrag_main/UniFrag/runUniFrag/zn_cr_cifs_noduplicated/cifs_heavy_elements"
    
    os.makedirs(target_dir, exist_ok=True)
    
    cif_paths = sorted(glob.glob(os.path.join(cif_dir, "*.cif")))
    print(f"Scanning {len(cif_paths)} parent CIFs in {cif_dir}...")
    
    heavy_elements = {'I', 'Si', 'Br', 'B', 'Se', 'As'}
    moved_count = 0
    moved_files = []
    
    for path in cif_paths:
        refcode = os.path.splitext(os.path.basename(path))[0].upper()
        try:
            reader = EntryReader(path)
            entry = reader[0]
            ccdc_mol = entry.molecule
            
            elements_in_cif = {atom.atomic_symbol for atom in ccdc_mol.atoms}
            matched_heavy = elements_in_cif.intersection(heavy_elements)
            
            if matched_heavy:
                # Move this file
                dest_path = os.path.join(target_dir, os.path.basename(path))
                shutil.move(path, dest_path)
                moved_count += 1
                moved_files.append((refcode, sorted(list(matched_heavy))))
        except Exception as e:
            print(f"Error reading parent {refcode}: {e}")
            
    print(f"\nPurification Complete!")
    print(f"Successfully moved {moved_count} CIF files containing heavy/semi-metal elements to:")
    print(f"  {target_dir}")
    print(f"Remaining CIF files in parent directory: {len(glob.glob(os.path.join(cif_dir, '*.cif')))}")
    
    print("\nMoved structures detail:")
    print("| REFCODE | Elements Found |")
    print("| :--- | :--- |")
    for ref, elems in sorted(moved_files):
        print(f"| **{ref}** | {', '.join(elems)} |")

if __name__ == "__main__":
    main()
