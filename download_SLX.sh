#! /bin/sh
#
# download_SLX.sh
# Copyright (C) 2016 macint01 <macint01@jblab-srv001>
# Copyright (C) 2021 Thomas Bradley <thomas.bradley@cruk.cam.ac.uk>

# Distributed under terms of the MIT license.
#
if [ $# -ne 2 ]; then
	    echo $0: usage: ./download_SLX.sh \<SLX ID\> \<output_directory\>
	    exit 1
fi

SLX=$1
OUTPUT_DIR=$2

/home/bioinformatics/software/java/java21/bin/java -jar /home/bioinformatics/software/clarity-tools.jar -l $SLX -f "*1.fq.gz"  --download-dir $OUTPUT_DIR
/home/bioinformatics/software/java/java21/bin/java -jar /home/bioinformatics/software/clarity-tools.jar -l $SLX -f "*2.fq.gz" --download-dir $OUTPUT_DIR
# /home/bioinformatics/software/java/java21/bin/java -jar /home/bioinformatics/software/clarity-tools.jar -l $SLX -f "*.ba*" --download-dir $OUTPUT_DIR
# /home/bioinformatics/software/java/java21/bin/java -jar /home/bioinformatics/software/clarity-tools.jar -l $SLX -f "*.tar*" --download-dir $OUTPUT_DIR
/home/bioinformatics/software/java/java21/bin/java -jar /home/bioinformatics/software/clarity-tools.jar -l $SLX -f "*contents*" --download-dir $OUTPUT_DIR
/home/bioinformatics/software/java/java21/bin/java -jar /home/bioinformatics/software/clarity-tools.jar -l $SLX -f "*html*" --download-dir $OUTPUT_DIR
