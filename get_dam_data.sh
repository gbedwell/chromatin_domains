#!/bin/bash

#$ -S /bin/bash
#$ -N align_DamID
#$ -pe pvm 20
#$ -l m_mem_free=12G
#$ -wd /home/gjb30/tsa_seq_hmm
#$ -j n
#$ -m bea
#$ -M gregoryj_bedwell@dfci.harvard.edu

source "/home/gjb30/miniforge3/etc/profile.d/conda.sh"
conda activate intmap

f_index_prefix="/mnt/storage/labs/aengelman/T2T/genome/indexes/bowtie2_noY/hs1_noY"
m_index_prefix="/mnt/storage/labs/aengelman/T2T/genome/indexes/bowtie2/hs1"
rscript="/mnt/storage/apps/R-4.4.1/bin/Rscript"
id=<access_key_id>
secret=<access_key_secret>

echo "Parsing LMNB1 DamID metadata..."

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

for target in "laminB1_DamID"; do

    echo "Cleaning "$target" metadata..."
    $rscript --verbose process_metadata.R "$target"

    echo "Downloading "$target" DamID-seq data..."
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
    
    echo "Cropping "$target" DamID-seq reads..."

    cd "$dvar"
    for input_file in *_replicate_*.fastq.gz; do
        filename=$(basename "$input_file" .fastq.gz)
        cutadapt \
            -j 20 \
            -e 0.1 \
            -m 10 \
            -g CGCGGCCGAGGA \
            -o ${filename}_trimmed.fastq.gz \
            $input_file
    done

    echo "Aligning "$target" DamID-seq reads..."

    for input_file in *_replicate_*_trimmed.fastq.gz; do
        filename=$(basename "$input_file" .fastq.gz)

        if [[ "$filename" == *"K562"* ]]; then
            index_prefix="$f_index_prefix"
        else
            index_prefix="$m_index_prefix"
        fi

        bowtie2 --very-sensitive -x "$index_prefix" -U "$input_file" -p 20 --no-unal |\
        samtools view -b -h -@ 20 - |\
        samtools sort -@ 20 -o "$filename"_rmdup.bam -
        samtools index -@ 20 -b "$filename"_rmdup.bam "$filename"_rmdup.bai

        rm "$filename".fastq.gz

    done
    cd ..
    echo
done

echo
echo "Done!"

conda deactivate