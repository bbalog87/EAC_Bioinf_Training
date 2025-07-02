#!/bin/bash

# Default values
DEFAULT_CPU_COUNT=8

# Function to display usage information
usage() {
    cat << EOM
Usage: $0 -i <input_reads_folder> -o <output_folder> [-c <cpu_count>]
Options:
  -i, --input-folder   Path to the folder containing raw read files (required)
  -o, --output-folder  Path to the output folder for processed files (required)
  -c, --cpu-count      Number of CPUs to use (optional, default: $DEFAULT_CPU_COUNT)
  -h, --help           Display this help message
EOM
    exit 1
}

# Parse command-line parameters
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-folder)
            READS_DIR="$2"
            shift
            shift
            ;;
        -o|--output-folder)
            FASTP_OUT="$2"
            shift
            shift
            ;;
        -c|--cpu-count)
            CPU_COUNT="$2"
            shift
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if required parameters are provided
if [ -z "$READS_DIR" ] || [ -z "$FASTP_OUT" ]; then
    usage
fi

# Set default CPU count if not provided
if [ -z "$CPU_COUNT" ]; then
    CPU_COUNT=$DEFAULT_CPU_COUNT
fi

# Check if the input reads folder exists
if [ ! -d "$READS_DIR" ]; then
    echo "Error: Input reads folder '$READS_DIR' does not exist."
    exit 1
fi

# Create the output folder if it doesn't exist
mkdir -p "$FASTP_OUT"

# Process raw reads using fastp
process_reads() {
    local fwd_file="$1"
    local base="$(basename "${fwd_file}" _1.fastq.gz)"
    local rev_file="${fwd_file%_1.fastq.gz}_2.fastq.gz"

    # Check if the reverse read file exists
    if [ ! -f "$rev_file" ]; then
        echo "Error: $rev_file does not exist."
        exit 1
    fi

    # Create output directory for each sample
    local sample_output="$FASTP_OUT/$base"
    mkdir -p "$sample_output"

    # Run fastp
    fastp -i "$fwd_file" -o "$sample_output/${base}.fastp_1.fastq" \
          -I "$rev_file" -O "$sample_output/${base}.fastp_2.fastq" \
          -t 5 --cut_tail 5 -q 28 -p -g \
          -j "$sample_output/${base}_fastp_report.json" \
          --html "$sample_output/${base}_fastp_report.html" \
          -R "Fastp analysis for reads data" -w "$CPU_COUNT" -l 40

    # Check for errors
    if [ $? -ne 0 ]; then
        echo "Error: fastp failed for $fwd_file and $rev_file"
        exit 1
    fi
}

# Iterate through raw read files
for fwd_file in "$READS_DIR"/*_1.fastq.gz; do
    process_reads "$fwd_file"
done


echo "=========================="
START_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "$START_DATE_TIME : Fastp processing completed successfully for all samples."
echo "=========================="
