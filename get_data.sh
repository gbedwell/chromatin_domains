#!/bin/bash

#$ -S /bin/bash
#$ -N align_tsa
#$ -pe pvm 20
#$ -l m_mem_free=12G
#$ -wd /home/gjb30/tsa_seq_hmm
#$ -j n
#$ -m bea
#$ -M gregoryj_bedwell@dfci.harvard.edu

source "/home/gjb30/miniforge3/etc/profile.d/conda.sh"
conda activate intmap

f_index_prefix="/mnt/storage/labs/aengelman/T2T/genome_noY_bt2/chm13v2.0_noY"
m_index_prefix="/mnt/storage/labs/aengelman/T2T/genome_bt2/chm13v2.0"
rscript="/mnt/storage/apps/R-4.4.1/bin/Rscript"
id=<access_key_id>
secret=<access_key_secret>

echo "Parsing SON, LMNB1 metadata..."

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

for target in "son" "laminB1"; do

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

        bowtie2 --very-sensitive -x "$index_prefix" -U "$input_file" -p 20 --no-unal |\
        samtools view -b -h -@ 20 - |\
        samtools collate -O -@ 20 - |\
        samtools fixmate -@ 20 -m - - |\
        samtools sort -@ 20 - |\
        samtools markdup -r -@ 20 - "$filename"_rmdup.bam
        samtools index -@ 20 -b "$filename"_rmdup.bam "$filename"_rmdup.bai

        rm "$filename".fastq.gz

    done
    cd ..
    echo
done

echo
echo "Done!"

conda deactivate