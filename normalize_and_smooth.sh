#!/bin/bash

python2="/mnt/storage/apps/Python-2.7.18/bin/python"

for target in "son" "laminB1" "laminB1_DamID"; do

    sub_dir="${target}_data"

    for filename in "$sub_dir"/*input*; do 
        if [[ -f "$filename" && "$filename" == *".bam" ]]; then 
            type=$(basename "$filename" | cut -d '_' -f 1)
            rep=$(basename "$filename" | cut -d '_' -f 4)

            name="${type}_replicate_${rep}"
            input="${type}_input_replicate_${rep}_rmdup.bam"
            pulldown="${type}_pulldown_replicate_${rep}_rmdup.bam"

            if [[ "$filename" == *"K562"* ]]; then
                chrom_sizes="../chm13v2.0_noM_noY_chrom.sizes"
            else
                chrom_sizes="../chm13v2.0_noM_chrom.sizes"
            fi

            $python2 normalize_TSA-seq.py \
                -r 1000 \
                -w 20000 \
                -e "${sub_dir}/${pulldown}" \
                -c "${sub_dir}/${input}" \
                -g "$chrom_sizes" \
                --output "${sub_dir}/${name}_normalized" \
                --wig2bw /mnt/storage/apps/UCSC-tools/wigToBigWig
        fi
    done

    for filename in "$sub_dir"/*normalized.bw; do
        name=$(basename $filename)
        prefix=${name%_*}

        if [[ "$name" == *"K562"* ]]; then
            chrom_sizes="../chm13v2.0_noM_noY_chrom.sizes"
        else
            chrom_sizes="../chm13v2.0_noM_chrom.sizes"
        fi

        $python2 TSA_smooth.py \
            --bw "${filename}" \
            -g "$chrom_sizes" \
            -w 25000 \
            --smooth \
            -o "${sub_dir}" \
            -n "${prefix}_smooth"
    done

    types=()
    for cell_type in "${sub_dir}"/*smooth.wig; do
        types+=($(basename "$cell_type" | cut -d '_' -f 1))
    done

    unique_types=($(printf "%s\n" "${types[@]}" | sort -u))

    for name in "${unique_types[@]}"; do
        foi=("${sub_dir}"/"${name}"*smooth.wig)

        python2 combine_replicates.py \
            --wig1 "${foi[0]}" \
            --wig2 "${foi[1]}"  \
            -n "${sub_dir}/${name}_combined"
    done
done