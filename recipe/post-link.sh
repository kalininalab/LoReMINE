#!/bin/bash

# Locate the IPA executable inside the conda environment
IPA_EXEC="$PREFIX/bin/ipa"   # or ipa.py or the exact file name

# Only continue if the file exists
if [[ -f "$IPA_EXEC" ]]; then
    # Remove the --reason element in the array
    sed -i "/--reason/d" "$IPA_EXEC"
fi

exit 0