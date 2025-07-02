#!/bin/bash

# Usage information
usage() {
    echo -e "\033[1;34m[INFO::] Usage: $0 -i <input_reads_folder> -p <project_folder> -d <kraken2_db> -q <krakenuniq_db> -k <kaiju_db_dir> -c <centrifuge_db> -r <reference_genome> [-t <threads>]\033[0m"
    echo ""
    echo "Options:"
    echo "  -i, --input          Path to the folder containing renamed FASTQ files (required)."
    echo "  -p, --project        Path to the project folder where all results will be stored (required)."
    echo "  -d, --kraken2_db     Path to the Kraken2 database (required)."
    echo "  -q, --krakenuniq_db  Path to the KrakenUniq database (required)."
    echo "  -k, --kaiju_db_dir   Path to the Kaiju database directory (required)."
    echo "  -c, --centrifuge_db  Path to the Centrifuge database index (required)."
    echo "  -r, --reference      Path to the reference genome for alignment and variant calling (required)."
    echo "  -t, --threads        Number of threads to use (optional, default: 60)."
    echo "  -h, --help           Display this help message."
    exit 1
}

# Default values
THREADS=14

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_DIR="$2"; shift 2;;
        -p|--project) PROJECT_DIR="$2"; shift 2;;
        -d|--kraken2_db) KRAKEN_DB="$2"; shift 2;;
        -q|--krakenuniq_db) KRAKENUNIQ_DB="$2"; shift 2;;
        -k|--kaiju_db_dir) KAIJU_DB_DIR="$2"; shift 2;;
        -c|--centrifuge_db) CENTRIFUGE_DB="$2"; shift 2;;
        -r|--reference) REFERENCE="$2"; shift 2;;
        -t|--threads) THREADS="$2"; shift 2;;
        -h|--help) usage;;
        *) echo "Unknown option: $1"; usage;;
    esac
done

# Check if required arguments are provided
if [[ -z "$INPUT_DIR" || -z "$PROJECT_DIR" || -z "$KRAKEN_DB" || -z "$KRAKENUNIQ_DB" || -z "$KAIJU_DB_DIR" || -z "$CENTRIFUGE_DB" || -z "$REFERENCE" ]]; then
    echo -e "\033[1;31m[ERROR::] Missing required arguments.\033[0m"
    usage
fi

# Create main output directories
KRAKEN_OUT="$PROJECT_DIR/Kraken2"
KRAKENUNIQ_OUT="$PROJECT_DIR/KrakenUniq"
KAIJU_OUT="$PROJECT_DIR/Kaiju"
CENTRI_OUT="$PROJECT_DIR/Centrifuge"
ASSEMBLY_OUT="$PROJECT_DIR/Assembly"
VARIANTS_OUT="$PROJECT_DIR/Variants"
ALIGNMENT_OUT="$PROJECT_DIR/Alignment"
QUAST_OUT="$PROJECT_DIR/Quast"
mkdir -p "$KRAKEN_OUT" "$KRAKENUNIQ_OUT" "$KAIJU_OUT" "$CENTRI_OUT" "$ASSEMBLY_OUT" "$VARIANTS_OUT" "$ALIGNMENT_OUT" "$QUAST_OUT"

# Loop through all _R1.fastq files and process
for fwd_file in "$INPUT_DIR"/*_R1.fastq; do
    sample=$(basename "$fwd_file" | sed -r 's/_R1.fastq//')
    rev_file="$INPUT_DIR/${sample}_R2.fastq"

    echo -e "\033[1;34m[INFO::] Processing sample $sample...\033[0m"
    
    cat "$fwd_file" "$rev_file" > sample.merged.fq

    # Kraken2 classification
    echo -e "\033[1;33m[INFO::] Running Kraken2 for $sample...\033[0m"
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" sample.merged.fq \
            --unclassified-out "$KRAKEN_OUT/${sample}_unclassified.txt" \
            --classified-out "$KRAKEN_OUT/${sample}_classified.txt" \
            --output "$KRAKEN_OUT/${sample}_kraken_out" \
            --report "$KRAKEN_OUT/${sample}_kraken_report.txt" --use-names  --memory-mapping   

    rm sample.merged.fq

    # Extract non-host (non-human) reads based on Kraken2 results
    awk '$1 == "C" && $3 != "Homo sapiens (taxid 9606)" {print $2}' \
        "$KRAKEN_OUT/${sample}_kraken_out" > "$KRAKEN_OUT/${sample}_nonHost_read_ids.txt"




# 1. Filter non-host reads from R1 and R2 using SeqKit
echo -e "\033[1;33m[INFO::] Filtering non-host reads for $sample...\033[0m"
seqkit grep -f "$KRAKEN_OUT/${sample}_nonHost_read_ids.txt" -o "$KRAKEN_OUT/${sample}_nonHost_R1.fastq" -j "$THREADS" "$INPUT_DIR/${sample}_R1.fastq"
seqkit grep -f "$KRAKEN_OUT/${sample}_nonHost_read_ids.txt" -o "$KRAKEN_OUT/${sample}_nonHost_R2.fastq" -j "$THREADS" "$INPUT_DIR/${sample}_R2.fastq"

# 2. Ensure perfectly paired reads with SeqKit and output to correctly named filtered files
seqkit pair -1 "$KRAKEN_OUT/${sample}_nonHost_R1.fastq" \
            -2 "$KRAKEN_OUT/${sample}_nonHost_R2.fastq" \
            -O "$KRAKEN_OUT" -j "$THREADS" -u 


 # Check if the filtered files were created
if [[ ! -f "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" || ! -f "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" ]]; then
    echo -e "\033[1;31m[ERROR::] Paired non-host reads were not created for $sample.\033[0m"
    exit 1
fi




# 3. Alignment to reference using BWA MEM
bwa index "$REFERENCE/mpxv_ref.fna"






# Inside the loop, for each sample being processed
RG_INFO="@RG\tID:${sample}\tSM:${sample}\tLB:lib1\tPL:ILLUMINA\tPU:${sample}.flowcell.1101"

echo -e "\033[1;33m[INFO::] Aligning filtered reads to reference for $sample...\033[0m"
bwa mem -t "$THREADS" -R "$RG_INFO"  "$REFERENCE/mpxv_ref.fna" \
       "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" \
       "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" | \
    samtools view -b - | samtools sort -o "$ALIGNMENT_OUT/${sample}_sorted.bam"

# Verify if the sorted BAM file was created successfully
if [[ ! -f "$ALIGNMENT_OUT/${sample}_sorted.bam" ]]; then
    echo -e "\033[1;31m[ERROR::] BAM file generation failed for $sample.\033[0m"
    exit 1
fi


   # Extract paired reads
    samtools view -b -f 1 -F 12 "$ALIGNMENT_OUT/${sample}_sorted.bam" > "$ALIGNMENT_OUT/${sample}_properly_paired.bam"
    samtools fastq -1 "$ALIGNMENT_OUT/${sample}_paired_R1.fastq" \
                   -2 "$ALIGNMENT_OUT/${sample}_paired_R2.fastq" \
                   "$ALIGNMENT_OUT/${sample}_properly_paired.bam"
    
 
 # Mark duplicates with new Picard syntax
    picard MarkDuplicates -I "$ALIGNMENT_OUT/${sample}_properly_paired.bam" \
                          -O "$ALIGNMENT_OUT/${sample}_dedup.bam" \
                          -M "$ALIGNMENT_OUT/${sample}_metrics.txt" \
                          --REMOVE_DUPLICATES true




    # Generate consensus and SNP calling
    samtools mpileup -aa -A -Q 0 -d 10000 -f "$REFERENCE/mpxv_ref.fna" \
                     "$ALIGNMENT_OUT/${sample}_dedup.bam" | \
                     ivar consensus -p "$VARIANTS_OUT/${sample}_consensus.fa" -m 10 -n N -t 0.75
    samtools mpileup -aa -A -Q 20 -d 10000 -f "$REFERENCE/mpxv_ref.fna" \
                    "$ALIGNMENT_OUT/${sample}_dedup.bam" | \
                    ivar variants -p "$VARIANTS_OUT/${sample}_variants" -m 30 -q 20 -t 0.75 > "$VARIANTS_OUT/${sample}_mutations.vcf"




# Ensure perfectly paired reads with SeqKit and output to correctly named filtered files
seqkit pair -1 "$ALIGNMENT_OUT/${sample}_paired_R1.fastq" \
            -2 "$ALIGNMENT_OUT/${sample}_paired_R2.fastq" \
            -O "$ALIGNMENT_OUT" -j "$THREADS" -u  --force
mv "$ALIGNMENT_OUT/${sample}_paired_R1.paired.fastq" "$ALIGNMENT_OUT/${sample}_paired_R1.fastq"
mv "$ALIGNMENT_OUT/${sample}_paired_R2.paired.fastq" "$ALIGNMENT_OUT/${sample}_paired_R2.fastq" 




    # KrakenUniq classification
    echo -e "\033[1;33m[INFO::] Running KrakenUniq for $sample...\033[0m"
    krakenuniq --db "$KRAKENUNIQ_DB" --paired \
               --threads "$THREADS" \
               --classified-out "$KRAKENUNIQ_OUT/${sample}_classified.txt" \
               --output "$KRAKENUNIQ_OUT/${sample}_krakenuniq_out" \
               --report-file "$KRAKENUNIQ_OUT/${sample}_report.tsv" \
               "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" 


    # Kraken2 classification with paired reads
    echo -e "\033[1;33m[INFO::] Running Kraken2 for $sample...\033[0m"

  kraken2 --db "$KRAKEN_DB" --paired --threads "$THREADS" \
        --output "$KRAKEN_OUT/${sample}_paired_kraken_out" \
        --report "$KRAKEN_OUT/${sample}_paired_kraken_report.txt" \
        --use-names \
        "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq"



    # Bracken processing on Kraken2 and KrakenUniq output
    echo -e "\033[1;33m[INFO::] Running Bracken for $sample...\033[0m"
    bracken -d "$KRAKEN_DB" -i "$KRAKEN_OUT/${sample}_paired_kraken_report.txt" \
            -o "$KRAKEN_OUT/${sample}_paired_bracken_output.txt" \
            -w "$KRAKEN_OUT/${sample}_paired_bracken_report.txt" -l S -t "$THREADS"

    bracken -d "$KRAKEN_DB" -i "$KRAKENUNIQ_OUT/${sample}_report.tsv" \
            -o "$KRAKENUNIQ_OUT/${sample}_bracken_output.txt" \
            -w "$KRAKENUNIQ_OUT/${sample}_uniq_bracken_report.txt" -l S -t "$THREADS"

   
# Kaiju classification
echo -e "\033[1;33m[INFO::] Running Kaiju for $sample...\033[0m"

# Step 1: Run Kaiju for each sample
kaiju -t "${KAIJU_DB_DIR}/nodes.dmp" \
      -f "${KAIJU_DB_DIR}/kaiju_db_viruses.fmi" \
      -i "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" \
      -j "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" \
      -z "$THREADS" -o "${KAIJU_OUT}/${sample}_kaiju.out"

# Step 2: Add taxon names to the Kaiju output
kaiju-addTaxonNames -t "${KAIJU_DB_DIR}/nodes.dmp" \
                    -n "${KAIJU_DB_DIR}/names.dmp" \
                    -i "${KAIJU_OUT}/${sample}_kaiju.out" \
                    -o "${KAIJU_OUT}/${sample}_kaiju-names.out"

# Step 3: Filter and format the classified results
grep '^C' "${KAIJU_OUT}/${sample}_kaiju-names.out" | sed 's/;/\t/g' > "${KAIJU_OUT}/${sample}_classified-kaiju.tmp"

# Step 4: Add sample name as a prefix to each line and save to final output
awk -v base="$sample" '{OFS="\t"; print base, $0}' "${KAIJU_OUT}/${sample}_classified-kaiju.tmp" > "${KAIJU_OUT}/${sample}_classified-kaiju.out"

# Step 5: Clean up temporary file
rm "${KAIJU_OUT}/${sample}_classified-kaiju.tmp"

# Print completion message
echo -e "\033[1;32m[INFO::] Kaiju processing completed for $sample.\033[0m"

# Centrifuge classification
echo -e "\033[1;33m[INFO::] Running Centrifuge for $sample...\033[0m"

if [ -f "$KRAKEN_OUT/${sample}.single.fastq" ]; then
    centrifuge -x "$CENTRIFUGE_DB" \
               -1 "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" \
               -2 "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" \
               -U "$KRAKEN_OUT/${sample}.single.fastq" \
               -p "$THREADS" \
               -S "${CENTRI_OUT}/${sample}_centrifuge_out" \
               --report-file "${CENTRI_OUT}/${sample}_centrifuge_report.txt"
else
    centrifuge -x "$CENTRIFUGE_DB" \
               -1 "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" \
               -2 "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" \
               -p "$THREADS" \
               -S "${CENTRI_OUT}/${sample}_centrifuge_out" \
               --report-file "${CENTRI_OUT}/${sample}_centrifuge_report.txt"
fi

               

    # De novo assembly with SPAdes
    echo -e "\033[1;33m[INFO::] Running SPAdes assembly for $sample...\033[0m"

cat $KRAKEN_OUT/${sample}_nonHost_R1.unpaired.fastq  $KRAKEN_OUT/${sample}_nonHost_R2.unpaired.fastq >$KRAKEN_OUT/${sample}.single.fastq
rm $KRAKEN_OUT/${sample}_nonHost_R?.unpaired.fastq 


  spades.py --isolate \
              -1 "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" \
              -2 "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq" \
              -s $KRAKEN_OUT/${sample}.single.fastq \
              --only-assembler -t "$THREADS" -o "$ASSEMBLY_OUT/$sample"
    
    # Assembly stats
    assembly-stats "$ASSEMBLY_OUT/$sample/contigs.fasta" > "$ASSEMBLY_OUT/all_asm_stats.txt"
	
	
	
	  # KrakenUniq classification for assembly 
    krakenuniq --db "$KRAKENUNIQ_DB"  \
               --threads "$THREADS" \
               --classified-out "$ASSEMBLY_OUT/$sample/${sample}_contigs_krakenUniq_classified.txt" \
               --output "$ASSEMBLY_OUT/$sample/${sample}_contigs_krakenuniq_out" \
               --report-file "$ASSEMBLY_OUT/$sample/${sample}_contigs_krakenUniq_report.tsv" \
               "$ASSEMBLY_OUT/$sample/contigs.fasta"
	
	  bracken -d "$KRAKEN_DB" -i "$ASSEMBLY_OUT/$sample/${sample}_contigs_krakenUniq_report.tsv" \
            -o "$ASSEMBLY_OUT/$sample/${sample}_contigs_brackenUniq_output.txt" \
            -w "$ASSEMBLY_OUT/$sample/${sample}_contigs_brackenUniq_report.txt" -l S -t "$THREADS"
			
			
	
	  kraken2 --db "$KRAKEN_DB" --threads "$THREADS" \
        --output "$ASSEMBLY_OUT/$sample/${sample}_contigs_kraken_out" \
        --report "$ASSEMBLY_OUT/$sample/${sample}_contigs_kraken_report.txt" \
        --use-names \
         "$ASSEMBLY_OUT/$sample/contigs.fasta"
		
	
	bracken -d "$KRAKEN_DB" -i "$ASSEMBLY_OUT/$sample/${sample}_contigs_kraken_report.txt" \
            -o "$ASSEMBLY_OUT/$sample/${sample}_contigs_bracken_output.txt" \
            -w "$ASSEMBLY_OUT/${sample}_contigs_bracken_report.txt" -l S -t "$THREADS"
		
	
	## Runing centrfige on contig assembly for each sample
	 centrifuge -x "$CENTRIFUGE_DB" -f \
                -U "$ASSEMBLY_OUT/$sample/contigs.fasta" \
               -p "$THREADS" \
               -S "$ASSEMBLY_OUT/$sample/${sample}_contigs_centrifuge_out" \
               --report-file "$ASSEMBLY_OUT/$sample/${sample}_contigs_centrifuge_report.txt"
	
	
	# Step 1: Run Kaiju on assembly copntoig fpore each sample
kaiju -t "${KAIJU_DB_DIR}/nodes.dmp" \
      -f "${KAIJU_DB_DIR}/kaiju_db_viruses.fmi" \
      -i "$ASSEMBLY_OUT/$sample/contigs.fasta" \
      -z "$THREADS" -o "$ASSEMBLY_OUT/$sample/${sample}_contigs_kaiju.out"

# Step 2: Add taxon names to the Kaiju output
kaiju-addTaxonNames -t "${KAIJU_DB_DIR}/nodes.dmp" \
                    -n "${KAIJU_DB_DIR}/names.dmp" \
                    -i "$ASSEMBLY_OUT/$sample/${sample}_contigs_kaiju.out" \
                    -o "$ASSEMBLY_OUT/$sample/${sample}_contigs_kaiju-names.out"
	
	
	
	
	
	
	
	
	
	
	
	
	
	

    # QUAST Assessment
   
 echo -e "\033[1;33m[INFO::] Running QUAST Asssessment for $sample...\033[0m"

# Index the contigs.fasta for alignment
bwa index "$ASSEMBLY_OUT/$sample/contigs.fasta"

# Run QUAST for assembly quality assessment
quast -o "$QUAST_OUT/${sample}" \
      -r "$REFERENCE/mpxv_ref.fna" -t "$THREADS" \
      --circos -g "$REFERENCE/mpxv_ref.gff" \
      "$ASSEMBLY_OUT/$sample/contigs.fasta" \
      --ref-bam "$ALIGNMENT_OUT/${sample}_sorted.bam" \
      -1 "$KRAKEN_OUT/${sample}_nonHost_R1.paired.fastq" \
      -2 "$KRAKEN_OUT/${sample}_nonHost_R2.paired.fastq"


		  
   # conda activate mpxv_env
 
    echo -e "\033[1;32m[INFO::] Finished processing $sample.\033[0m"
done

# Generate a comprehensive Kaiju summary table after processing all samples

cat ${KAIJU_OUT}/S*_classified-kaiju.out | awk '$2 == "C"' >${PROJECT_DIR}/allSamples_kaiju_results.txt
cat "$ASSEMBLY_OUT/$sample/S*_kaiju-names.out" | awk '$2 == "C"' >${PROJECT_DIR}/allSamples_kaiju_contigs_results.txt

# Check if Kaiju output files exist
if ls "${KAIJU_OUT}"/S*_kaiju-names.out 1> /dev/null 2>&1; then
    echo -e "[INFO::] Generating Kaiju summary table for all samples...\n"
	
     kaiju2table -c 200 -u \
            -o ${PROJECT_DIR}/allSamples_kaiju.table"\
			-t "${KAIJU_DB_DIR}/nodes.dmp" \
			-n "${KAIJU_DB_DIR}/names.dmp" \
			-r species -e ${KAIJU_OUT}/S*_kaiju-names.out
			
			
	    kaiju2table -c 200 -u \
            -o ${PROJECT_DIR}/allSamples_contigs_kaiju.table" \
			-t "${KAIJU_DB_DIR}/nodes.dmp" \
			-n "${KAIJU_DB_DIR}/names.dmp" \
			-r species -e "$ASSEMBLY_OUT/$sample/S*_contigs_kaiju-names.out"
					
		
			
    echo -e "[INFO::] Kaiju summary table saved as ${PROJECT_DIR}/allSamples_kaiju_summary.table \n"
	echo -e "[INFO::] Kaiju summary table saved as ${PROJECT_DIR}/allSamples_contigs_kaiju_summary.table \n"
	
    awk 'NR==1 {$1="sample"; print} NR>1 {match($1, /S[0-9]+/); \
    print substr($1, RSTART, RLENGTH), $2, $3, $4, $5}' ${PROJECT_DIR}/allSamples_kaiju.table" > ${PROJECT_DIR}/allSamples_kaiju_summary.table"
	
	awk 'NR==1 {$1="sample"; print} NR>1 {match($1, /S[0-9]+/); \
    print substr($1, RSTART, RLENGTH), $2, $3, $4, $5}' ${PROJECT_DIR}/allSamples_contigs_kaiju_summary.table > ${PROJECT_DIR}/allSamples_contigs_kaiju_summary.table"

else
    echo -e "[ERROR::] No Kaiju output files found in ${KAIJU_OUT}.\n"
fi

echo -e "\033[1;32m[INFO::] Pipeline completed for all samples!\033[0m"



 









