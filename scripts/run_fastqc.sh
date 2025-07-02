#!/usr/bin/env bash

# Author: Julien Nguinkal
# Date: [2023-01-20]
# Description: This is a wrapper to perform initial QC on raw reads using the FastQC command line tool

# Exit if any command exits with a non-zero status
set -e
# Exit if any variable is unset
set -u

# Exit if any command in a pipeline exits with a non-zero status
set -o pipefail

#### Color on terminal output

BBlue='\033[1;34m' 
RED='\033[0;31m'
NOCOLOR='\033[0m'

### Check user inputs
CPUS=8
OUT=./FastQC_output
FASTQC_BINARY=fastqc

usage() {
    echo "Usage: $0 -d <directory> -t <threads> -o <output_directory> -f <fastqc_binary>"
    echo " -d <directory> : specify the directory containing paired reads"
    echo " -t <threads> : specify the number of threads to use for fastqc (default: 8)"
    echo " -o <output_directory> : specify the output directory for fastqc results (default: ./FastQC_output)"
    echo " -f <fastqc_binary> : specify the path to fastqc binary (default: fastqc)"
    exit 1
}



### Check user inputs

while getopts ":d:t:o:f:h" opt; do

    case $opt in
        d)
            READS_DIR="$OPTARG"
            ;;
        t)
            CPUS="$OPTARG"
            ;;
        o)
            OUT="$OPTARG"
            ;;
        f)
            FASTQC_BINARY="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
		h)
            usage
            ;;
    esac
done


if [[ "${@:$OPTIND}" =~ "-h" || "${@:$OPTIND}" =~ "--help" ]]; then
    usage
fi


if [ $# -eq 0 ]
then
    echo "Use -d <directory> to specify the directory containing paired reads"
	 echo "Use -h or -help to print extensisve options"
    exit 1
fi


## Print start of the process and date time
echo "==========================" 
START_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e  "${BBlue}$START_DATE_TIME : RUNNING ...fastqc for all samples...${NOCOLOR}" 
 echo "==========================" 

# Check if the READS_DIR directory exists and contains fastq files
if [ ! -d "$READS_DIR" ] || [ -z "$(ls -A $READS_DIR)" ]; then
    echo "Error: $READS_DIR directory does not exist or is empty."
    exit 1
fi


# Check if the output directory exists
if [ -d "$OUT" ]; then
    echo "Warning: Output directory already exists. Files may be overwritten."
else
    mkdir -p $OUT
fi

# Check if fastqc is installed and in the PATH or specified by user
if [ -z "$FASTQC_BINARY" ]; then
    if ! command -v fastqc > /dev/null; then
        echo "Error: fastqc is not installed or not in the PATH."
        exit 1
    fi
else
    if ! command -v "$FASTQC_BINARY" > /dev/null; then
        echo "Error: fastqc binary not found at $FASTQC_BINARY"
        exit 1
    fi
fi

# Run fastqc to check reads quality - Visual inspection needed.
for file in $READS_DIR/*_1.fastq.gz; do
    base=$(basename "$file" _1.fastq.gz)
    #base=${file%_1.fastq.gz} ## sample name
    fwd=$file ## forward reads
    rev=${file%_1.fastq.gz}_2.fastq.gz ## reverse reads
    # Check if the reverse read file exists
    if [[ ! -f "$rev" ]]; then
        echo "Error: $rev does not exist."
        exit 1
    fi
    mkdir -p $OUT/$base/ ## output dir for each sample
    $FASTQC_BINARY $fwd $rev -t $CPUS -o $OUT/$base/
    # check for errors
    if [ $? -ne 0 ]; then
        echo "Error: fastqc failed for $fwd and $rev"
        exit 1
    fi
done

echo "==========================" 
START_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e  "${BBlue}$START_DATE_TIME : fastqc processing completed successfully for all samples. Have fun!${NOCOLOR}" 
echo "==========================" 


