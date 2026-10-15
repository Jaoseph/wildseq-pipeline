# a whitelist of 10X cell barcodes present in the sample is prepared from the transcriptome library
# this requires umi-tools to be installed
# enrich_pcr.sh

#!/bin/bash


UNMAPPED_TXT=$1
UNMAPPED_COUNT=$2
SLX_ID=$3
SAMPLE_INDEX=$4
ENRICH_INDEX=$5
TABLE_FULL=$6
TABLE_FULL_FILTERED=$7

WHITELIST_STDIN="data/raw_data/${SLX_ID}/${SAMPLE_INDEX}*R1*.fastq.gz"
EXTRACT_STDIN="data/raw_data/${SLX_ID}/${ENRICH_INDEX}*R1*.fastq.gz"
EXTRACT_READ2IN="data/raw_data/${SLX_ID}/${ENRICH_INDEX}*R2*.fastq.gz"

# define dynamic paths for intermediate files
WHITELIST_OUT="results/${SLX_ID}/${SAMPLE_INDEX}/enrich_pcr/whitelist_${SAMPLE_INDEX}_mRNA_All.txt"
EXTRACTED_R2="results/${SLX_ID}/${SAMPLE_INDEX}/enrich_pcr/${ENRICH_INDEX}_PCR_All_R2_extracted.fastq.gz"
BARCODE_FASTA="results/${SLX_ID}/${SAMPLE_INDEX}/enrich_pcr/${ENRICH_INDEX}_PCR_All_R2_extracted_barcode_sequences_only.fasta"
PCR_BOWTIE="results/${SLX_ID}/${SAMPLE_INDEX}/enrich_pcr/${ENRICH_INDEX}_PCR.bowtie"

#comment missing.
mkdir -p results/${SLX_ID}/${SAMPLE_INDEX}/enrich_pcr

umi_tools whitelist --stdin $WHITELIST_STDIN --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNN \
--log2stderr > $WHITELIST_OUT

umi_tools extract --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin $EXTRACT_STDIN \
--read2-stdout --read2-in $EXTRACT_READ2IN \
--stdout $EXTRACTED_R2 \
--whitelist $WHITELIST_OUT

# the part of the read matching our WILDseq barcode pattern is extracted
zcat $EXTRACTED_R2 | \
egrep -A 2 -B 1 [ATGC]{12}GCTACCTGGTCGTTCTATCG[ATGC]{12} | \
sed '/^--$/d' | sed -n '1~4s/^@/>/p;2~4p'| \
sed 's/\(.*\)\([ATCG]\{12\}GCTACCTGGTCGTTCTATCG[ATCG]\{12\}\)\(.*\)/\2/g' > $BARCODE_FASTA

# the barcode sequences are aligned to an index of all possible barcode sequences based on the origin
#twist oligo pool from which the vector library was generated
software/bowtie-1.3.0/bowtie -f -v 2 --un $UNMAPPED_TXT \
-x data/theoretical_barcodes2/WILDseq_index2 $BARCODE_FASTA > $PCR_BOWTIE

sort $UNMAPPED_TXT | uniq -c | sort -nrk 1 > $UNMAPPED_COUNT


# the 10X cell barcode, UMI and WILDseq barcode alignment information is collated for each read and
#duplicates collapsed - this is as has always been done previously
awk -F'["\t"_ ]' '{print $2 "\t" $3 "\t" $6}' $PCR_BOWTIE | sort | uniq > $TABLE_FULL

# decided might be better to further collapse any reads with the same CBC and UMI but a different WILDseq
#barcode. This version will only keep WILDseq barcode with the highest count for each unique CBC/UMI
#combination. Where ties exist all barcodes are listed
awk -F'["\t"_ ]' '{print $2 "\t" $3 "\t" $6}' $PCR_BOWTIE | sort | uniq -c | \
sed 's/^[ \t]*//' | awk -F'["\t" ]' '{print $1 "\t" $2 "\t" $3 "\t" $4}' | \
awk '{k=$2FS$3; if (a[k]<$1){a[k]=$1; data[k]=$0} else if (a[k]==$1) {a[k]=$1; data[k]=data[k]","$4}} END {for (i in a) print data[i]}' \
> $TABLE_FULL_FILTERED

