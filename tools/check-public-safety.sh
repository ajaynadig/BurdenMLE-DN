#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

fail() {
  echo "PUBLIC-SAFETY CHECK FAILED: $*" >&2
  exit 1
}

# Manuscript inputs and outputs must never enter the public repository.
for forbidden_dir in analysis/inputs analysis/outputs analysis/logs; do
  [[ ! -e "${forbidden_dir}" ]] ||
    fail "restricted/generated directory exists: ${forbidden_dir}"
done

# Disallow data-bearing formats except the explicitly synthetic package file
# and the metadata-only input manifest.
while IFS= read -r path; do
  case "${path}" in
    ./inst/extdata/synthetic_example.csv|./analysis/input_manifest.tsv)
      ;;
    *)
      fail "unexpected data-like file: ${path}"
      ;;
  esac
done < <(
  find . -type f \( \
    -iname '*.csv' -o -iname '*.tsv' -o -iname '*.txt' -o \
    -iname '*.vcf' -o -iname '*.vcf.gz' -o -iname '*.bcf' -o \
    -iname '*.bed' -o -iname '*.bed.gz' -o -iname '*.bgz' -o \
    -iname '*.rds' -o -iname '*.rdata' -o -iname '*.rda' \
  \) -print
)

# Absolute user paths are a portability and disclosure risk.
if rg -n '/Users/|/home/|~/Mirror/' \
    --glob '!tools/check-public-safety.sh' . >/dev/null; then
  rg -n '/Users/|/home/|~/Mirror/' \
    --glob '!tools/check-public-safety.sh' .
  fail "absolute user path found"
fi

# High-signal filename terms catch accidental additions before a commit.
if find . -type f \
    \( -iname '*variant*table*' -o -iname '*count*table*' -o \
       -iname '*pedigree*' -o -iname '*participant*' \) \
    ! -path './analysis/scripts/*' -print | grep -q .; then
  find . -type f \
    \( -iname '*variant*table*' -o -iname '*count*table*' -o \
       -iname '*pedigree*' -o -iname '*participant*' \) \
    ! -path './analysis/scripts/*' -print
  fail "suspicious data filename found"
fi

echo "Public-safety check passed: no manuscript count/variant data detected."
