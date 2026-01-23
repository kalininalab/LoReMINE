#!/bin/bash

# Locate the IPA executable inside the conda environment
IPA_EXEC="$PREFIX/bin/ipa"   # or ipa.py or the exact file name

# Only continue if the file exists
if [[ -f "$IPA_EXEC" ]]; then
    # Remove the --reason element in the array
    sed -i "/--reason/d" "$IPA_EXEC"
fi

# Locate download_files.py anywhere inside the opt folder (handles version changes like dfast_qc-1.0.7)
DFAST_FILE=$(find "$PREFIX/opt" -name "download_files.py" -path "*/dqc/download_files.py" | head -n 1)

# Only continue if the file was found
if [[ -n "$DFAST_FILE" && -f "$DFAST_FILE" ]]; then
    # Replace os.makedirs(out_dir) with os.makedirs(out_dir, exist_ok=True)
    sed -i "s/os.makedirs(out_dir)/os.makedirs(out_dir, exist_ok=True)/" "$DFAST_FILE"
fi

exit 0