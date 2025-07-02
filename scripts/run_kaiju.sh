#!/bin/bash

# Usage function
usage() {
  echo "Usage: $(basename "$0") -r <reads_directory> -g <host_genome> [-t <threads>] [-h]"
  echo "  -r: Path to the reads directory."
  echo "  -g: Path to the host reference genome."
  echo "  -t: Number of threads to use. Default is 8."
  echo "  -h: Display this help message."
  exit 1
}

# Default values
threads=8

# Parse command-line options
while getopts "r:g:t:h" opt; do
  case $opt in
    r) reads_dir="$OPTARG";;
    g) reference_genome="$OPTARG";;
    t) threads="$OPTARG";;
    h) usage;;
    \?) echo "Invalid option: -$OPTARG"; usage;;
  esac
done

# Check if required options are provided
if [ -z "$reads_dir" ] || [ -z "$reference_genome" ]; then
  echo "Error: Both -r and -g options are required."
  usage
fi

# Verify input values
if [ ! -d "$reads_dir" ] || [ ! -f "$reference_genome" ]; then
  echo "Error: Input data not found. Make sure reads directory and reference genome are valid."
  usage
fi

echo "Reads Directory: $reads_dir"
echo "Reference Genome: $reference_genome"

# Set paths
SCRIPT="/home/nguinkal/Scripts"
WORKING_DIR=$(pwd)
KAIJU_DIR="$WORKING_DIR/kaiju_dir"
FASTP_OUT_DIR="$WORKING_DIR/fastp_out" # Directory for trimmed reads

echo "Working Directory: $WORKING_DIR"

# Create working directory
mkdir -p "$KAIJU_DIR"
cd "$KAIJU_DIR" || exit

# Extract the filename from the provided path
reference_genome_filename=$(basename "$reference_genome")
reference_genome_filename_nonHost="${reference_genome_filename}_nonHost"

# Build bowtie2 index
mkdir -p host.genome
cp "$reference_genome" host.genome

echo "Building Bowtie2 Index from: $reference_genome"

#bowtie2-build "$reference_genome" \host.genome/"$reference_genome_filename" --threads "$threads"

# Iterate through trimmed read files in the fastp_out directory
for fwd_file in "$FASTP_OUT_DIR"/*_fastp_1.fastq.gz; do
    if [[ -f "$fwd_file" ]]; then
        # Determine base name
        base=$(basename "$fwd_file" _fastp_1.fastq.gz)
        
        # Create output folder
        mkdir -p "$base"
      #  cd "$base" || exit
        
        # Determine paths for forward and reverse reads
        fwd_read="$fwd_file"
        rev_read="$FASTP_OUT_DIR/${base}_fastp_2.fastq.gz"
        
        echo "Mapping $base: Forward: $fwd_read, Reverse: $rev_read"
        
        # Mapping with bowtie2 using trimmed reads
       # bowtie2 -1 "$fwd_read" -2 "$rev_read" \
          #      -S "$base.sam" --un-conc "nonHost/${base}_reads_unmapped.fastq" \
              #  --threads "$threads" -x "$WORKING_DIR/host.genome/$reference_genome_filename"

        # Move back to the working directory
        
    fi
done

echo "Metagenomics analysis pipeline completed!"
