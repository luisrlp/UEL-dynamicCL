#!/bin/bash
# Add the MUGRAD and JFLUX vector fields (add_vector_field.py) to every ODB of a study folder.
#
# Usage (from test_in_abaqus/):
#   ./add_vector_field_all.sh                  # all ODBs in SA_stage2/, one at a time
#   ./add_vector_field_all.sh SA_stage2 4      # 4 ODBs at a time
#   ./add_vector_field_all.sh SA_stage3        # another study folder
#   ABAQUS=abq2024 ./add_vector_field_all.sh   # other Abaqus command
#
# Each ODB is opened for writing and saved in place. Frames that already have the fields are
# skipped by add_vector_field.py, so the script can be rerun. Jobs still running (.lck file
# present) are skipped. The output of each ODB goes to <run folder>/add_vector_field.log, and
# the list of results to <study folder>/add_vector_field_all.log.

SA_DIR="${1:-SA_stage2}"
NPAR="${2:-1}"
SCRIPT="$(cd "$(dirname "$0")" && pwd)/add_vector_field.py"
ABAQUS="${ABAQUS:-abaqus}"

[ -d "$SA_DIR" ] || { echo "Folder $SA_DIR not found"; exit 1; }
[ -f "$SCRIPT" ] || { echo "$SCRIPT not found"; exit 1; }

process() {
    local odb="$1" dir job
    dir=$(dirname "$odb")
    job=$(basename "$odb" .odb)
    if [ -e "$dir/$job.lck" ]; then
        echo "SKIPPED (running)  $odb"
        return 0
    fi
    if (cd "$dir" && "$ABAQUS" python "$SCRIPT" "$job" < /dev/null > add_vector_field.log 2>&1) \
        && grep -q "Script completed successfully" "$dir/add_vector_field.log"; then
        echo "OK                 $odb"
    else
        echo "FAILED             $odb (see $dir/add_vector_field.log)"
    fi
}
export -f process
export SCRIPT ABAQUS

LOG="$SA_DIR/add_vector_field_all.log"
N=$(find "$SA_DIR" -name '*.odb' ! -name 'upgraded_*' | wc -l)
echo "Adding vector fields to $N ODB files in $SA_DIR ($NPAR at a time)..."
find "$SA_DIR" -name '*.odb' ! -name 'upgraded_*' | sort \
    | xargs -r -P "$NPAR" -I{} bash -c 'process "$1"' _ {} | tee "$LOG"
echo "Done: $(grep -c '^OK' "$LOG") OK, $(grep -c '^FAILED' "$LOG") failed, $(grep -c '^SKIPPED' "$LOG") skipped (list in $LOG)"
