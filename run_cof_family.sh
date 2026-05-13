#!/bin/bash
# Run UniFrag COF fragmentation for every CIF in a folder.
# Usage: ./run_cof_family.sh <folder> [radius] [auto|monomer|dimer|both]
# Default layer mode is both, so monomer and dimer files are saved separately.
# Example: ./run_cof_family.sh test_on_cof_zn_pc_series 4.0

set -u

ROOT="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-/Users/omert/miniconda3/bin/python}"
SCRIPT="$ROOT/fragmentation_oop.py"
FOLDER_ARG="${1:-}"
RADIUS="${2:-4.0}"
COF_LAYER="${3:-both}"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

if [ -z "$FOLDER_ARG" ]; then
    echo "Usage: $0 <folder> [radius] [auto|monomer|dimer|both]"
    echo "Example: $0 test_on_cof_zn_pc_series 4.0"
    exit 1
fi

if [[ "$FOLDER_ARG" = /* ]]; then
    FAMILY_DIR="$FOLDER_ARG"
else
    FAMILY_DIR="$ROOT/$FOLDER_ARG"
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

SCRIPT_HELP="$($PYTHON $SCRIPT -h 2>&1 || true)"
SUPPORTS_COF_LAYER=0
if echo "$SCRIPT_HELP" | grep -q -- "--cof-layer"; then
    SUPPORTS_COF_LAYER=1
fi

case "$COF_LAYER" in
    auto|monomer|dimer|both) ;;
    *)
        echo -e "${RED}Invalid COF layer:${NC} $COF_LAYER"
        echo "Use auto, monomer, dimer, or both."
        exit 1
        ;;
esac

if [ "$SUPPORTS_COF_LAYER" -ne 1 ] && [ "$COF_LAYER" != "auto" ]; then
    echo -e "${YELLOW}Warning:${NC} current fragmentation_oop.py does not support --cof-layer."
    echo "  Falling back to auto layer mode for compatibility."
    COF_LAYER="auto"
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

run_variant() {
    local cif="$1"
    local base="$2"
    local layer="$3"
    local suffix=""

    if [ "$layer" != "auto" ]; then
        suffix="_${layer}"
    fi

    local out_norm="$FAMILY_DIR/${base}_frag_cof${suffix}.xyz"
    local out_min="$FAMILY_DIR/${base}_frag_cof${suffix}_min.xyz"

    local normal_output
    if [ "$layer" = "auto" ] || [ "$SUPPORTS_COF_LAYER" -ne 1 ]; then
        normal_output=$("$PYTHON" "$SCRIPT" "$cif" --kind cof --radius "$RADIUS" --output "$out_norm" 2>&1)
    else
        normal_output=$("$PYTHON" "$SCRIPT" "$cif" --kind cof --radius "$RADIUS" --cof-layer "$layer" --output "$out_norm" 2>&1)
    fi
    local normal_status=$?
    echo "$normal_output" | grep -E "COF Path|Path [A-Z]|Final size|Saved|Error|Warning|Traceback" || true

    local min_output
    if [ "$layer" = "auto" ] || [ "$SUPPORTS_COF_LAYER" -ne 1 ]; then
        min_output=$("$PYTHON" "$SCRIPT" "$cif" --kind cof --radius "$RADIUS" --minimize --output "$out_min" 2>&1)
    else
        min_output=$("$PYTHON" "$SCRIPT" "$cif" --kind cof --radius "$RADIUS" --cof-layer "$layer" --minimize --output "$out_min" 2>&1)
    fi
    local min_status=$?
    echo "$min_output" | grep -E "COF Path|Path [A-Z]|Final size|Saved|Error|Warning|Traceback" || true

    if [ $normal_status -ne 0 ] || [ $min_status -ne 0 ]; then
        echo -e "  ${RED}FAILED${NC} ${layer}"
        if [ $normal_status -ne 0 ]; then
            failures+="$base $layer normal failed\n"
            echo "$normal_output"
        fi
        if [ $min_status -ne 0 ]; then
            failures+="$base $layer minimum failed\n"
            echo "$min_output"
        fi
        return 1
    fi

    local norm_info min_info norm_atoms norm_formula min_atoms min_formula
    norm_info="$(formula_for_xyz "$out_norm")"
    min_info="$(formula_for_xyz "$out_min")"
    norm_atoms="$(echo "$norm_info" | cut -f1)"
    norm_formula="$(echo "$norm_info" | cut -f2-)"
    min_atoms="$(echo "$min_info" | cut -f1)"
    min_formula="$(echo "$min_info" | cut -f2-)"

    summary+="$base\t$layer\t$norm_atoms\t$norm_formula\t$min_atoms\t$min_formula\n"
    echo -e "  ${GREEN}OK${NC} ${layer}: normal=$norm_atoms min=$min_atoms"
    return 0
}

run_one() {
    local cif="$1"
    local base
    base="$(basename "${cif%.*}")"

    echo -e "${YELLOW}${base}${NC}"

    local failed=0
    if [ "$COF_LAYER" = "both" ]; then
        run_variant "$cif" "$base" monomer || failed=1
        run_variant "$cif" "$base" dimer || failed=1
    else
        run_variant "$cif" "$base" "$COF_LAYER" || failed=1
    fi

    if [ $failed -ne 0 ]; then
        ((fail++))
    else
        ((pass++))
    fi
}

echo "============================================================"
echo "  UniFrag COF Family Test"
echo "============================================================"
echo "  Folder: $FAMILY_DIR"
echo "  Radius: $RADIUS"
echo "  COF layer: $COF_LAYER"
echo "  CIF files: ${#cifs[@]}"
echo ""

for cif in "${cifs[@]}"; do
    run_one "$cif"
    echo ""
done

echo "============================================================"
echo -e "name\tlayer\tnormal_atoms\tnormal_formula\tmin_atoms\tmin_formula"
printf "%b" "$summary"
echo "============================================================"
echo -e "Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"
if [ $fail -ne 0 ]; then
    echo "Failures:"
    printf "%b" "$failures"
    exit 1
fi
exit 0
