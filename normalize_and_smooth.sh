python2="/mnt/storage/apps/Python-2.7.18/bin/python"
sub_dir="son_data"
for filename in "$sub_dir"/*input*; do 
    if [[ -f "$filename" && "$filename" == *".bam" ]]; then 
        type=$(basename "$filename" | cut -d '_' -f 1)
        rep=$(basename "$filename" | rev | cut -d '_' -f 2 | rev )
        name="${type}_replicate_${rep}"
        input="${type}_input_replicate_${rep}_rmdup.bam"
        pulldown="${type}_pulldown_replicate_${rep}_rmdup.bam"

        $python2 normalize_TSA-seq2.0.py \
            -r 20000 \
            -w 20000 \
            -e "${sub_dir}/${pulldown}" \
            -c "${sub_dir}/${input}" \
            -g ../chm13v2.0_noM_chrom.sizes \
            --output "${sub_dir}/${name}_normalized" \
            --wig2bw /mnt/storage/apps/UCSC-tools/wigToBigWig
    fi
done

for filename in "$sub_dir"/*normalized.wig; do
    name=$(basename $filename)
    prefix=${name%_*}

    $python2 TSA_smooth_hanningFor20kbNonsliding_TSA2.0.py \
        --wig "${filename}" \
        -w 20000 \
        -aggwin 200000 \
        --smooth \
        -n1 "${sub_dir}/${prefix}_sm_20kb" \
        -n2 "${sub_dir}/${prefix}_sm_20kb_avg_200kb" \
        -g ../chm13v2.0_noM_chrom.sizes
done

sub_dir="laminB1_data"
for filename in "$sub_dir"/*input*; do
    if [[ -f "$filename" && "$filename" == *".bam" ]]; then
        type=$(basename "$filename" | cut -d '_' -f 1)
        rep=$(basename "$filename" | rev | cut -d '_' -f 2 | rev )
        name="${type}_replicate_${rep}"
        input="${type}_input_replicate_${rep}_rmdup.bam"
        pulldown="${type}_pulldown_replicate_${rep}_rmdup.bam"

        $python2 normalize_TSA-seq2.0.py \
            -r 20000 \
            -w 20000 \
            -e "${sub_dir}/${pulldown}" \
            -c "${sub_dir}/${input}" \
            -g ../chm13v2.0_noM_chrom.sizes \
            --output "${sub_dir}/${name}_normalized" \
            --wig2bw /mnt/storage/apps/UCSC-tools/wigToBigWig
    fi
done

for filename in "$sub_dir"/*normalized.wig; do
    name=$(basename $filename)
    prefix=${name%_*}

    $python2 TSA_smooth_hanningFor20kbNonsliding_TSA2.0.py \
        --wig "${filename}" \
        -w 20000 \
        -aggwin 200000 \
        --smooth \
        -n1 "${sub_dir}/${prefix}_sm_20kb" \
        -n2 "${sub_dir}/${prefix}_sm_20kb_avg_200kb" \
        -g ../chm13v2.0_noM_chrom.sizes
done
