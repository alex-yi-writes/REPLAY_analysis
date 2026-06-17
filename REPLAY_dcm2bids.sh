#!/bin/bash

# paths
raw_data_dir="/Users/yyi/Desktop/REPLAY/sourcedata"
output_dir="/Users/yyi/Desktop/REPLAY/rawdata"
# config_file="/Users/yyi/Desktop/REPLAY/codes/REPLAY_dcm2bids.json"
config_file="/Users/yyi/Desktop/REPLAY/codes/REPLAY_postupdate_dcm2bids.json"
log_dir="/Users/yyi/Desktop/REPLAY/logs"

mkdir -p "${log_dir}"

# --- Subject list: 102 to xxx -------------------------------
subject_ids=()
for i in $(seq 147); do
    subject_ids+=("${i}")
done

# --- Conversion loop ----------------------------------------
for subject in "${subject_ids[@]}"; do

    echo "·̩̩̥͙＊*•̩̩͙✩•̩̩͙*˚　　˚*•̩̩͙✩•̩̩͙*˚＊·̩̩̥͙"
    echo "Processing subject ${subject}"
    echo "·̩̩̥͙＊*•̩̩͙✩•̩̩͙*˚　　˚*•̩̩͙✩•̩̩͙*˚＊·̩̩̥͙"

    subject_raw_dir="${raw_data_dir}/${subject}"

    if [[ ! -d "${subject_raw_dir}" ]]; then
        echo "WARNING: sourcedata directory not found for subject ${subject}, skipping."
        continue
    fi

    # Locate session folder one level below the subject directory
    # (sourcedata/102/<session_dir>/NNN_series_name/)
    session_dir=$(ls -d "${subject_raw_dir}"/*/  2>/dev/null | head -1)

    if [[ -z "${session_dir}" ]]; then
        echo "WARNING: no session directory found for subject ${subject}, skipping."
        continue
    fi

    echo "  source: ${session_dir}"

    dcm2bids \
        -d "${session_dir}" \
        -p "${subject}" \
        -c "${config_file}" \
        -o "${output_dir}" \
        # --forceDcm2niix \
        # --clobber \
        2>&1 | tee "${log_dir}/sub-${subject}_dcm2bids.log"

    exit_code=${PIPESTATUS[0]}
    if [[ ${exit_code} -ne 0 ]]; then
        echo "ERROR: dcm2bids failed for subject ${subject} (exit code ${exit_code})"
        echo "  see ${log_dir}/sub-${subject}_dcm2bids.log"
    fi

    # Gzip any uncompressed NIfTIs (belt-and-braces; dcm2niix should already produce .nii.gz)
    find "${output_dir}/sub-${subject}" -name "*.nii" -exec gzip -f {} \;

    echo "finished subject ${subject}"
    echo ""

done

echo "all subjects done"
