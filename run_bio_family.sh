#!/bin/bash
# Run UniFrag bio-molecule sliding-window fragmentation for every PDB in a folder.
# Usage:   ./run_bio_family.sh <folder> [window_size] [stride] [ph]
# Example: ./run_bio_family.sh test_on_bio_mol 5 1 7.0
# Pass --no-pdbfixer as 5th arg to skip PDBFixer:
#        ./run_bio_family.sh test_on_bio_mol 5 1 7.0 --no-pdbfixer

set -u

ROOT="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-/Users/omert/miniconda3/bin/python}"
SCRIPT="$ROOT/fragmentation_oop.py"
FOLDER_ARG="${1:-}"
WINDOW="${2:-5}"
STRIDE="${3:-1}"
PH="${4:-7.0}"
EXTRA_FLAG="${5:-}"          # pass --no-pdbfixer here if desired

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# Usage guard
# ---------------------------------------------------------------------------
if [ -z "$FOLDER_ARG" ]; then
    echo "Usage: $0 <folder> [window_size] [stride] [ph] [--no-pdbfixer]"
    echo "Example: $0 test_on_bio_mol 5 1 7.0"
    exit 1
fi

if [[ "$FOLDER_ARG" = /* ]]; then
    FAMILY_DIR="$FOLDER_ARG"
else
    FAMILY_DIR="$ROOT/$FOLDER_ARG"
fi

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
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
pdbs=("$FAMILY_DIR"/*.pdb)
if [ ${#pdbs[@]} -eq 0 ]; then
    echo -e "${RED}No PDB files found in:${NC} $FAMILY_DIR"
    exit 1
fi

# ---------------------------------------------------------------------------
# Helper: read window count and atom stats from the summary CSV
# ---------------------------------------------------------------------------
csv_summary() {
    "$PYTHON" - "$1" <<'PYSUM'
import csv, sys
from pathlib import Path
p = Path(sys.argv[1])
if not p.exists():
    print("0\t0\t0\t0")
    sys.exit(0)
rows = list(csv.DictReader(p.open()))
if not rows:
    print("0\t0\t0\t0")
    sys.exit(0)
n_win   = len(rows)
totals  = [int(r["n_total"])  for r in rows]
heavies = [int(r["n_heavy"])  for r in rows]
print(f"{n_win}\t{min(totals)}\t{max(totals)}\t{round(sum(totals)/n_win,1)}")
PYSUM
}

# ---------------------------------------------------------------------------
# Per-PDB runner
# ---------------------------------------------------------------------------
pass=0
fail=0
summary=""
failures=""

run_one() {
    local pdb="$1"
    local base
    base="$(basename "${pdb%.*}")"
    local out_dir="$FAMILY_DIR/bio_fragments_${base}"

    echo -e "${YELLOW}${base}${NC}"

    local run_output
    run_output=$("$PYTHON" "$SCRIPT" "$pdb" \
        --kind bio \
        --window-size "$WINDOW" \
        --stride      "$STRIDE" \
        --ph          "$PH" \
        --output-dir  "$out_dir" \
        ${EXTRA_FLAG} \
        2>&1)
    local run_status=$?

    # Print key lines (suppress noisy RDKit UFF warnings)
    echo "$run_output" | grep -v "UFFTYPER\|WARNING: could not find" \
        | grep -E "Window|Summary CSV|Total windows|Error|Traceback" || true

    if [ $run_status -ne 0 ]; then
        echo -e "  ${RED}FAILED${NC}"
        failures+="$base failed\n"
        echo "$run_output"
        ((fail++))
        return
    fi

    local csv_file="$out_dir/${base}_windows.csv"
    local stats
    stats="$(csv_summary "$csv_file")"
    local n_win min_atoms max_atoms avg_atoms
    n_win="$(echo "$stats"    | cut -f1)"
    min_atoms="$(echo "$stats" | cut -f2)"
    max_atoms="$(echo "$stats" | cut -f3)"
    avg_atoms="$(echo "$stats" | cut -f4)"

    summary+="${base}\t${n_win} windows\tatoms: ${min_atoms}-${max_atoms} (avg ${avg_atoms})\n"
    echo -e "  ${GREEN}OK${NC}: windows=${n_win}  atoms=${min_atoms}-${max_atoms} (avg ${avg_atoms})"
    ((pass++))
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  UniFrag Bio-Molecule Family Test"
echo "============================================================"
echo "  Folder     : $FAMILY_DIR"
echo "  Window size: $WINDOW residues"
echo "  Stride     : $STRIDE residues"
echo "  pH         : $PH"
echo "  PDBFixer   : ${EXTRA_FLAG:-(enabled)}"
echo "  PDB files  : ${#pdbs[@]}"
echo ""

for pdb in "${pdbs[@]}"; do
    run_one "$pdb"
    echo ""
done

echo "============================================================"
echo -e "name\twindows\tatom_range"
printf "%b" "$summary"
echo "============================================================"
echo -e "Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"

if [ $fail -ne 0 ]; then
    echo "Failures:"
    printf "%b" "$failures"
    exit 1
fi
exit 0
