mkdir -p data_copy

for target in "son" "laminB1" "laminB1_DamID"; do

    sub_dir="${target}_data"

    mkdir -p data_copy/${sub_dir}
    
    echo "Copying ${target} data..."
    cd ${sub_dir}
    cp *.wig "../data_copy/${sub_dir}"
    cp *.bw "../data_copy/${sub_dir}"
    cp *.bed "../data_copy/${sub_dir}"
    cd ..

done