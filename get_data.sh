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
f_index_prefix="/mnt/storage/labs/aengelman/T2T/genome_noY_bwa/chm13v2.0_noY.fa.gz"
m_index_prefix="/mnt/storage/labs/aengelman/T2T/genome_bwa/chm13v2.0.fa.gz"
rscript="/mnt/storage/apps/Downloads/R-4.2.1/bin/Rscript"
id="<access_key_id>"
secret="<access_key_secret>"

echo "Parsing SON, LMNB1,and LMNB1 DamID metadata..."

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

# H1-hESC replicates removed: 4DNFIK4DQEP6, 4DNFIEE9OW9A
# HCT116 replicates removed: 4DNFIJIBLGAV, 4DNFI5FS93PB, 4DNFIMWU62GW, 4DNFIX78LMRJ 
# HFFc6 replicates removed: 4DNFIRHCKHG1, 4DNFIDJF8TIS
# K562 replicates removed: 4DNFI6QZFIBH, 4DNFIECFZRE2, 4DNFIBFA1K6U, 4DNFIFDNGQ47
cat laminB1_DamID_metadata_2024-02-27-16h-17m.tsv |\
    grep -v ^# |\
    grep -v 4DNFIJIBLGAV |\
    grep -v 4DNFI5FS93PB |\
    grep -v 4DNFIMWU62GW |\
    grep -v 4DNFIX78LMRJ |\
    grep -v 4DNFIK4DQEP6 |\
    grep -v 4DNFIEE9OW9A |\
    grep -v 4DNFIRHCKHG1 |\
    grep -v 4DNFIDJF8TIS |\
    grep -v 4DNFI6QZFIBH |\
    grep -v 4DNFIECFZRE2 |\
    grep -v 4DNFIBFA1K6U |\
    grep -v 4DNFIFDNGQ47 |\
    grep fastq > laminB1_DamID_cleaned.tsv  

for target in "son" "laminB1" "laminB1_DamID"; do

    echo "Cleaning "$target" metadata..."
    $rscript --verbose process_metadata.R "$target"

    echo "Downloading "$target" TSA-seq data..."
    dvar="${target}_data"

    mkdir -p "$dvar"

    xargs -n 1 -I {} -P 8 sh -c "cd \"$dvar\"/ && curl -O -L --user \"$id\":\"$secret\" {}" < "$target"_tsa_download_links.txt

    #modified from https://unix.stackexchange.com/a/646449
    while read -r accession sample; do
        for name in "$dvar"/"$accession".*; do
            [ -e "$name" ] || continue
            experiment="$dvar"/"$sample.fastq.gz"
            mv -- "$name" "$experiment"
        done
    done < "$target"_tsa_file_match.txt

    echo "Aligning "$target" TSA-seq reads..."

    cd "$dvar"
    for input_file in *_replicate_*.fastq.gz; do
        filename=$(basename "$input_file" .fastq.gz)
        if [[ "$filename" == *"K562"* ]]; then
            index_prefix="$f_index_prefix"
        else
            index_prefix="$m_index_prefix"
        fi
        $bwa_path mem -t 8 $index_prefix "$input_file" > "$filename".sam
        $samtools_path view -bS -@ 8 "$filename".sam > "$filename".bam
        rm "$filename.sam"
        $samtools_path sort -@ 8 "$filename".bam -o "$filename"_sort.bam
        $samtools_path markdup -r "$filename"_sort.bam "$filename"_rmdup.bam
        $samtools_path index "$filename"_rmdup.bam
        rm "$filename".fastq.gz; rm "$filename".bam; rm "$filename"_sort.bam
    done
    cd ..
    echo
done

echo
echo "Done!"