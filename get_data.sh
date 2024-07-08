#!/bin/bash

#$ -S /bin/bash
#$ -N align_tsa
#$ -pe pvm 20
#$ -l m_mem_free=12G
#$ -wd /home/gjb30/tsa_seq_hmm
#$ -j n
#$ -m bea
#$ -M gregoryj_bedwell@dfci.harvard.edu

bt2_path="/mnt/storage/apps/bowtie2-2.5.2/bowtie2"
samtools_path="/mnt/storage/apps/samtools/1.9/bin/samtools"
f_index_prefix="/mnt/storage/labs/aengelman/T2T/genome_noY_bt2/chm13v2.0_noY"
m_index_prefix="/mnt/storage/labs/aengelman/T2T/genome_bt2/chm13v2.0"
rscript="/mnt/storage/apps/Downloads/R-4.2.1/bin/Rscript"
id=<access_key_id>
secret=<access_key_secret>

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

        if [[ "$target" == *"DamID" ]]; then
            # Global alignment with DamID data yielded very low alignment rates (~10-20%).
            # Switching to local yielded rates more comparable to global alignment of the TSA-seq data.
            $bt2_path --local --very-sensitive-local -x "$index_prefix" -U "$input_file" -p 20 --no-unal |\
            $samtools_path view -b -h -@ 20 - |\
            $samtools_path sort -@ 20 - -o "$filename"_rmdup.bam
            $samtools_path index -@ 20 "$filename"_rmdup.bam
        else
            $bt2_path --very-sensitive -x "$index_prefix" -U "$input_file" -p 20 --no-unal |\
            $samtools_path view -b -h -@ 20 - |\
            $samtools_path collate -O -@ 20 - |\
            $samtools_path fixmate -@ 20 -m - - |\
            $samtools_path sort -@ 20 - |\
            $samtools_path markdup -r -@ 20 - "$filename"_rmdup.bam
            $samtools_path index -@ 20 "$filename"_rmdup.bam
        fi

        rm "$filename".fastq.gz

    done
    cd ..
    echo
done

echo
echo "Done!"