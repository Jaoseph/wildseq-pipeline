#!/bin/bash

INPUT_BAM=$1
SLX_ID=$2
SAMPLE_INDEX=$3
BOWTIE_INDEX=$4
UNMAPPED_TXT=$5
UNMAPPED_COUNT=$6
CB_UB_TXT=$7
BOWTIE_OUTPUT=$8

WILDSEQ_BAM="results/${SLX_ID}/${SAMPLE_INDEX}/barcode_extraction/WILDseq.bam"
BAMTOFASTQ_DIR="results/${SLX_ID}/${SAMPLE_INDEX}/barcode_extraction/bamtofastq"
FILTERED_FA="results/${SLX_ID}/${SAMPLE_INDEX}/barcode_extraction/bamtofastq_R2_filtered.fa"
CUT_FA="results/${SLX_ID}/${SAMPLE_INDEX}/barcode_extraction/bamtofastq_R2_filtered_cut.fa"


# Step 1: Extract reads mapping to the specific region
samtools view -b $INPUT_BAM \
EF1a_tdTomato_WILDseq:5398-5441 > $WILDSEQ_BAM

# Step 2: Convert BAM to FASTQ
software/cellranger-9.0.1/lib/bin/bamtofastq $WILDSEQ_BAM $BAMTOFASTQ_DIR

# Step 3: Filter reads containing the full barcode sequence
FILE=$(find results/${SLX_ID}/${SAMPLE_INDEX}/barcode_extraction/bamtofastq -name "*R2*" -exec realpath {} \;)
zcat ${FILE} | egrep -A 2 -B 1 [ATGC]{12}GCTACCTGGTCGTTCTATCG[ATGC]{12} | sed '/^--$/d'| \
sed -n '1~4s/^@/>/p;2~4p' > $FILTERED_FA

# Step 4: Trim flanking regions
cat $FILTERED_FA | sed 's/\(.*\)\([ATCG]\{12\}GCTACCTGGTCGTTCTATCG[ATCG]\{12\}\)\(.*\)/\2/g' > $CUT_FA

# Step 5: Align to theoretical barcode index
software/bowtie-1.3.0/bowtie -f -v 2 -x $BOWTIE_INDEX \
--un $UNMAPPED_TXT $CUT_FA > $BOWTIE_OUTPUT

# Step 6: Count unmapped barcodes
sort $UNMAPPED_TXT | uniq -c | sort -nrk 1 > $UNMAPPED_COUNT

# Step 7: Extract cell barcode and UMI info
samtools view $WILDSEQ_BAM | \
sed -r -n 's/.*(LH00417[:A-Z0-9]*).*(CB:Z:[ACTG]*).*(UB:Z:[ACTG]*).*/\1 \2 \3/p' > $CB_UB_TXT