#!/bin/bash

#$ -S /bin/bash
#$ -N mb_random
#$ -pe pvm 50
#$ -l m_mem_free=8G
#$ -wd /home/gjb30/xInt/large_random
#$ -j n
#$ -m bea
#$ -M gregoryj_bedwell@dfci.harvard.edu

/mnt/storage/apps/Downloads/R-4.4.1/bin/Rscript --verbose make_random.R &&
cat random_*.fa.gz > all.fa.gz

filename="all.fa.gz"
min_size_gb=1

bytes_to_gb() {
    local bytes=$1
    echo "scale=2; $bytes / 1024 / 1024 / 1024" | bc
}

if [ -f "$filename" ]; then
    file_size=$(stat -c %s "$filename")

    file_size_gb=$(bytes_to_gb "$file_size")

    if (( $(echo "$file_size_gb > $min_size_gb" | bc -l) )); then
        rm random_*.fa.gz
    else
        echo "Something went wrong."
    fi
else
    echo "$filename does not exist."
fi