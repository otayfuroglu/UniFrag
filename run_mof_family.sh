#!/bin/bash
# Run UniFrag MOF fragmentation for every CIF in a folder.
# Usage: ./run_mof_family.sh <folder> [radius]
# Example: ./run_mof_family.sh test_on_irmof_series 4.0

set -u

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-/Users/omert/miniconda3/bin/python}"
SCRIPT="/Users/omert/Desktop/UniFrag_main/UniFrag/fragmentation_oop.py"
FOLDER_ARG="${1:-}"
RADIUS="${2:-4.0}"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

if [ -z "$FOLDER_ARG" ]; then
    echo "Usage: $0 <folder> [radius]"
    echo "Example: $0 test_on_irmof_series 4.0"
    exit 1
fi

if [[ "$FOLDER_ARG" = /* ]]; then
    FAMILY_DIR="$FOLDER_ARG"
else
    FAMILY_DIR="$WORKDIR/$FOLDER_ARG"
fi

if [ ! -f "$PYTHON" ]; then
    echo -e "${RED}Missing Python:${NC} $PYTHON"
    exit 1
fi
if [ ! -f "$SCRIPT" ]; then
    echo -e "${RED}Missing script:${NC} $SCRIPT"
    exit 1
fi
if [ ! -d "$FAMILY_DIR" ]; then
    echo -e "${RED}Missing folder:${NC} $FAMILY_DIR"
    exit 1
fi

shopt -s nullglob
cifs=("$FAMILY_DIR"/*.cif)
if [ ${#cifs[@]} -eq 0 ]; then
    echo -e "${RED}No CIF files found in:${NC} $FAMILY_DIR"
    exit 1
fi

pass=0
fail=0
summary=""
csv_summary=""
failures=""

formula_for_xyz() {
    "$PYTHON" - "$1" <<'PYSUM'
from collections import Counter
from pathlib import Path
import sys
path = Path(sys.argv[1])
if not path.exists():
    print("MISSING\tMISSING")
    sys.exit(0)
species = []
for line in path.read_text().splitlines()[2:]:
    parts = line.split()
    if len(parts) >= 4:
        species.append(parts[0])
c = Counter(species)
formula = " ".join(f"{el}{c[el]}" for el in sorted(c))
print(f"{len(species)}\t{formula}")
PYSUM
}

run_one() {
    local cif="$1"
    local base
    base="$(basename "${cif%.*}")"
    local out_norm="$FAMILY_DIR/${base}_frag_mof.xyz"
    local out_min="$FAMILY_DIR/${base}_frag_mof_min.xyz"

    echo -e "${YELLOW}${base}${NC}"

    rm -f "$out_norm" "$out_min"

    local normal_output
    normal_output=$("$PYTHON" "$SCRIPT" "$cif" --kind mof --radius "$RADIUS" --output "$out_norm" 2>&1)
    local normal_status=$?
    echo "$normal_output" | grep -E "MOF Path|Path [A-Z]|Final size|Saved|Error|Warning|Traceback" || true

    if [ $normal_status -ne 0 ]; then
        echo -e "  ${RED}FAILED${NC}"
        failures+="$base normal failed\n"
        echo "$normal_output"
        ((fail++))
        return
    fi

    local norm_info min_info norm_atoms norm_formula min_atoms min_formula
    norm_info="$(formula_for_xyz "$out_norm")"
    norm_atoms="$(echo "$norm_info" | cut -f1)"
    norm_formula="$(echo "$norm_info" | cut -f2-)"

    # The python script automatically names the min output by removing .xyz and appending _min.xyz
    local expected_min_out="${out_norm%.xyz}_min.xyz"
    if [ -f "$expected_min_out" ]; then
        min_info="$(formula_for_xyz "$expected_min_out")"
        min_atoms="$(echo "$min_info" | cut -f1)"
        min_formula="$(echo "$min_info" | cut -f2-)"
    else
        min_atoms="N/A"
        min_formula="N/A"
    fi

    summary+="$base\t$norm_atoms\t$norm_formula\t$min_atoms\t$min_formula\n"
    csv_summary+="${base}.cif,$norm_atoms,$norm_formula,$min_atoms,$min_formula\n"
    echo -e "  ${GREEN}OK${NC}: normal=$norm_atoms min=$min_atoms"
    ((pass++))
}

echo "============================================================"
echo "  UniFrag MOF Family Test"
echo "============================================================"
echo "  Folder: $FAMILY_DIR"
echo "  Radius: $RADIUS"
echo "  CIF files: ${#cifs[@]}"
echo ""

for cif in "${cifs[@]}"; do
    run_one "$cif"
    echo ""
done

echo "============================================================"
echo -e "name\tnormal_atoms\tnormal_formula\tmin_atoms\tmin_formula"
printf "%b" "$summary"
echo "============================================================"

csv_file="$FAMILY_DIR/fragmentation_summary.csv"
echo "cif_file,normal_atoms,normal_formula,min_atoms,min_formula" > "$csv_file"
printf "%b" "$csv_summary" >> "$csv_file"
echo -e "CSV summary saved to: ${GREEN}$csv_file${NC}"
echo "============================================================"

echo -e "Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"
if [ $fail -ne 0 ]; then
    echo "Failures:"
    printf "%b" "$failures"
    exit 1
fi
exit 0
