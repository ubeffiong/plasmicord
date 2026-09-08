#!/usr/bin/env bash
# Shared helpers sourced by every stage script.
log()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2; }
die()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required command '$1' not found in PATH. Activate the conda env?"; }
have() { command -v "$1" >/dev/null 2>&1; }

# Iterate real rows of metadata.tsv. Prints TAB-separated fields keyed by header.
# We only strictly need: isolate_id, assembly_path. chromosomal_cluster/date/
# location are optional and passed through to later stages via the file itself.
read_isolates() {
    local sheet="$1"
    awk -F'\t' '
        NR==1 {
            for (i=1;i<=NF;i++) col[$i]=i
            if (!("isolate_id" in col))  { print "ERR:isolate_id" > "/dev/stderr"; exit 3 }
            if (!("assembly_path" in col)){ print "ERR:assembly_path" > "/dev/stderr"; exit 3 }
            next
        }
        /^[[:space:]]*#/ {next}
        /^[[:space:]]*$/ {next}
        {
            print $(col["isolate_id"]) "\t" $(col["assembly_path"])
        }
    ' "$sheet"
}
