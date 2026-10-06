#!/usr/bin/env bash
# Optional independent re-execution of the pinned official Linux binary.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo 'Usage: bash validate_ldblockshow.sh INPUT.vcf.gz OUTPUT_DIRECTORY' >&2
  exit 2
fi
input_path=$(realpath "$1")
output_path=$(realpath -m "$2")
mkdir -p "$output_path"
source_commit=1fc8cecbabf25416fa90de9c91a4e10656463738
archive_path="$output_path/LDBlockShow_pinned.tar.gz"
curl --fail --location --retry 2 "https://api.github.com/repos/hewm2008/LDBlockShow/tarball/$source_commit" --output "$archive_path"
mkdir -p "$output_path/source"
tar -xzf "$archive_path" -C "$output_path/source" --strip-components=1
chmod u+x "$output_path/source/bin/LDBlockShow" "$output_path/source/bin/ShowLDSVG"
sha256sum "$input_path" "$output_path/source/bin/LDBlockShow" > "$output_path/input_and_binary.sha256"
"$output_path/source/bin/LDBlockShow" -InVCF "$input_path" \
  -OutPut "$output_path/synthetic" -Region chr1:1000000:1100000 \
  -SeleVar 3 -BlockType 5 -MAF 0.01 -Miss 1 -Het 1 \
  > "$output_path/run.log" 2>&1
printf '%s\n' "$source_commit" > "$output_path/source_commit.txt"
echo "External phase-statistic outputs saved in $output_path"
