#!/bin/bash

# get flow cell name
tmp=${PWD##*/}
tmp=${tmp#*_*_*_}
FLOWCELL=${tmp%_*}

DUPLEX_CALLS_DIR="dorado_sup_duplex_bact_model"
mkdir -p $DUPLEX_CALLS_DIR

INPUT_DIR="./fastq_pass"
POD5_DIR="./pod5/"

### basecall with dorado for 5 kHz

# stop doradod service to use standalone dorado
sudo systemctl stop doradod.service 

# collect barcode dirs
dirs=$(find $INPUT_DIR -mindepth 1 -maxdepth 1 -type d -name '*barcode*')

for D in $dirs ; do
    echo ""
    echo "###############################################"
    
    barcode=${D##*/}
    echo "barcode is $barcode"

    ### get read ids from fastqs per barcode and duplex call those
    eval "$(conda shell.bash hook)"
    conda activate seqkit

    echo "collecting reads for $D"

    FASTQS=$(find $INPUT_DIR/$barcode -mindepth 1 -maxdepth 1 -type f -name '*.fastq.gz')
    
    for FASTQ in $FASTQS ; do
        seqkit fx2tab $FASTQ | awk -v OFS='\t' '{array[$1]=1} END {for (readID in array) print readID}' > $FASTQ\_reads.txt
    done

    cat $INPUT_DIR/$barcode/*.fastq.gz_reads.txt > $INPUT_DIR/$barcode/all_reads.txt
    rm $INPUT_DIR/$barcode/*.fastq.gz_reads.txt

    eval "$(conda shell.bash hook)"
    conda activate base

    mkdir -p $DUPLEX_CALLS_DIR/$barcode

    echo "basecalling data for $D"

    ~/ont/dorado/bin/dorado \
    duplex \
    /home/vc/ont/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@2023-09-22_bacterial-methylation \
    --read-ids $INPUT_DIR/$barcode/all_reads.txt \
    $POD5_DIR \
    -r \
    --min-qscore 10 \
    -x "cuda:0" \
    > ./$DUPLEX_CALLS_DIR/$barcode/duplex_calls.bam    

    ### split simplex and duplex
    duplex_only_dir="./fastq_pass_duplex_only/$barcode"
    duplex_simplex_dir="./fastq_pass_duplex_with_simplex/$barcode"
    
    mkdir -p $duplex_only_dir
    mkdir -p $duplex_simplex_dir

    eval "$(conda shell.bash hook)"
    conda activate samtools

    # duplex only
    samtools view -d dx:1 ./$DUPLEX_CALLS_DIR/$barcode/duplex_calls.bam > $duplex_only_dir/duplex_reads.bam
    samtools bam2fq $duplex_only_dir/duplex_reads.bam > $duplex_only_dir/$FLOWCELL\_pass_$barcode\_duplex_reads.fastq
    rm $duplex_only_dir/duplex_reads.bam

    # zip  fastqs
    pigz -f $duplex_only_dir/*.fastq

    # duplex and simplex
    # copy duplex over
    cp $duplex_only_dir/$FLOWCELL\_pass_$barcode\_duplex_reads.fastq.gz $duplex_simplex_dir

    samtools view -d dx:0 ./$DUPLEX_CALLS_DIR/$barcode/duplex_calls.bam > $duplex_simplex_dir/simplex_reads.bam
    samtools bam2fq $duplex_simplex_dir/simplex_reads.bam > $duplex_simplex_dir/$FLOWCELL\_pass_$barcode\_simplex_reads.fastq

    #cleanu up
    rm $duplex_simplex_dir/simplex_reads.bam

    # zip  fastqs
    pigz -f $duplex_simplex_dir/*.fastq

    eval "$(conda shell.bash hook)"
    conda activate base

done

# remove temp duplex calls dir
rm -rf ./$DUPLEX_CALLS_DIR

#restart doradod service
sudo systemctl start doradod.service 
