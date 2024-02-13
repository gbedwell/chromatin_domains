#!/bin/bash

#$ -S /bin/bash
#$ -N align_tsa
#$ -pe pvm 8
#$ -l m_mem_free=24000M
#$ -wd /mnt/storage/labs/aengelman/T2T/annotations/tsa_seq/
#$ -j n
#$ -m bea
#$ -M gregoryj_bedwell@dfci.harvard.edu

bwa_path="/mnt/storage/apps/BWA/bwa-0.7.17/bwa"
samtools_path="/mnt/storage/apps/samtools/1.9/bin/samtools"
index_prefix="/mnt/storage/labs/aengelman/T2T/genome_bwa/chm13v2.0.fa.gz"
rscript="/mnt/storage/apps/Downloads/R-4.2.1/bin/Rscript"
id="<access_key_id>"
secret="<access_key_id>"

echo "Parsing SON and LMNB1 TSA-seq metadata..."

#accession number 4DNFIUSJM8BY, labeled 'gDNA control', is K562 input control replicate 1
#grep fastq: only interested in fastq files
#grep -v 4DNFIAVMK1UI: remove one of the HFFc6 pulldown technical replicates (there are two; the choice was arbitrary)
#grep -v 4DNFIHMP7N14; remove one of the HCT116 pulldown technical replicates (arbitrary)
cat son_metadata_2024-01-08-09h-57m.tsv |\
    grep -v ^# |\
    grep fastq |\
    grep -v 4DNFIAVMK1UI |\
    grep -v 4DNFIHMP7N14 > son_cleaned.tsv

#grep -v 4DNFIICFO2J6: remove one of the HFFc6 input technical replicates (arbitrary)
#grep -v 4DNFI8OFH6ZV: remove one of the HFFc6 pulldown technical replicates (arbitrary)
cat laminB1_metadata_2024-01-08-09h-21m.tsv |\
    grep -v ^# |\
    grep fastq |\
    grep -v 4DNFIICFO2J6 |\
    grep -v 4DNFI8OFH6ZV > laminB1_cleaned.tsv

$rscript --verbose process_metadata.R

echo "Downloading SON TSA-seq data..."

mkdir -p son_data

xargs -n 1 -I {} -P 8 sh -c 'cd son_data/ && curl -O -L --user $id:$secret {}' < son_tsa_download_links.txt

#modified from https://unix.stackexchange.com/a/646449
while read -r accession sample; do
    for name in son_data/"$accession".*; do
        [ -e "$name" ] || continue
        experiment=son_data/"$sample.fastq.gz"
        mv -- "$name" "$experiment"
    done
done < son_tsa_file_match.txt

echo "Aligning SON TSA-seq reads..."

cd son_data
for input_file in *_replicate_*.fastq.gz; do
    filename=$(basename "$input_file" .fastq.gz)
    $bwa_path mem -t 8 $index_prefix "$input_file" > "$filename".sam
    $samtools_path view -bS -@ 8 "$filename".sam > "$filename".bam
    rm "$filename.sam"
    $samtools_path sort -@ 8 "$filename".bam -o "$filename"_sort.bam
    $samtools_path markdup -r "$filename"_sort.bam "$filename"_rmdup.bam
    $samtools_path index "$filename"_rmdup.bam
    rm "$filename".fastq.gz; rm "$filename".bam; rm "$filename"_sort.bam
done
cd ..

echo "Downloading LMNB1 TSA-seq data..."

mkdir -p laminB1_data

xargs -n 1 -I {} -P 8 sh -c 'cd laminB1_data/ && curl -O -L --user $id:$secret {}' < laminB1_tsa_download_links.txt

while read -r accession sample; do
    for name in laminB1_data/"$accession".*; do
        [ -e "$name" ] || continue
        experiment=laminB1_data/"$sample.fastq.gz"
        mv -- "$name" "$experiment"
    done
done < laminB1_tsa_file_match.txt

echo "Aligning LMNB1 TSA-seq reads..."

cd laminB1_data
for input_file in *_replicate_*.fastq.gz; do
    filename=$(basename "$input_file" .fastq.gz)
    $bwa_path mem -t 8 $index_prefix "$input_file" > "$filename".sam
    $samtools_path view -bS -@ 8 "$filename".sam > "$filename".bam
    rm "$filename.sam"
    $samtools_path sort -@ 8 "$filename".bam -o "$filename"_sort.bam
    $samtools_path markdup -r "$filename"_sort.bam "$filename"_rmdup.bam
    $samtools_path index "$filename"_rmdup.bam
    rm "$filename".fastq.gz; rm "$filename".bam; rm "$filename"_sort.bam
done
cd ..

echo "Done!"
