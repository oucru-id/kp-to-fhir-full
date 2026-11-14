#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 [OUTPUT_DIR] [CPUS]

Download the cgMLST.org Klebsiella pneumoniae/variicola/quasipneumoniae scheme
(2,358 targets, https://www.cgmlst.org/ncs/schema/2187931/) and convert it
into a chewBBACA schema with chewBBACA PrepExternalSchema.

  OUTPUT_DIR  default: db/cgmlst_schema (must not exist yet)
  CPUS        default: 4

Needs curl, unzip and chewBBACA.py 3.4.2 on PATH (bioconda::chewbbaca=3.4.2),
about 120 MB of download and 4 GB of disk. Then run the pipeline with
--cgmlst_schema OUTPUT_DIR. chewBBACA adds inferred alleles to the schema
during AlleleCall, so keep one copy per analysis series.
EOF
}

case "${1:-}" in -h|--help) usage; exit 0;; esac
out="${1:-db/cgmlst_schema}"
cpus="${2:-4}"
url="https://www.cgmlst.org/ncs/schema/2187931/alleles/"

if [ -e "$out" ]; then
    echo "ERROR: $out already exists" >&2
    exit 1
fi
for tool in curl unzip chewBBACA.py; do
    command -v "$tool" > /dev/null || { echo "ERROR: $tool not found on PATH" >&2; exit 1; }
done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
echo "Downloading $url"
curl -fL --retry 3 -o "$tmp/alleles.zip" "$url"
mkdir "$tmp/alleles"
unzip -q "$tmp/alleles.zip" -d "$tmp/alleles"
echo "$(find "$tmp/alleles" -name '*.fasta' | wc -l) allele files"

chewBBACA.py PrepExternalSchema -g "$tmp/alleles" -o "$out" --cpu "$cpus"
echo "chewBBACA schema written to $out"
