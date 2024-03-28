## Analyzing TSA-seq data with Multivariate HMMs

This repository contains the scripts necessary for identifying chromatin states from TSA-seq data with multivariate hidden Markov models. 

The data normalization steps are described in [Zhang, et. al.](10.1101/gr.266239.120). However, the normalized data were smoothed and binned according to [Wang, et. al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02253-3). TSA-seq and DamID-seq data were processed in the same way. The Python scripts used for these purposes were accordingly retrieved from the associated GitHub repositories ([Zhang](https://github.com/lgchang27/TSA-Seq-2020v2) and [Wang](https://github.com/zocean/Norma/tree/ab3153cb178933abfb03397a33de920d5d48b26a)).

The script <code>get_data.sh</code> downloads the required data from the [4DN data portal](https://data.4dnucleome.org/browse/?experimentset_type=replicate&type=ExperimentSetReplicate). A 4DN access key id and access key secret are required for this. The script then aligns the reads to the defined genome build using BWA MEM and performs further processing with samtools.

The script <code>normalize_and_smooth.sh</code> takes the sorted BAM files generated with <code>get_data.sh</code> and normalizes the reads, smooths them, and bins the smoothed data. Replicate data are then combined (averaged) to produce the values that go into the HMM.

The smoothed and binned data are then fed into <code>global_hmm.R</code> (in <code>funs.R</code>). The workhorse functionality of this function is importing the data directly from provided paths and automatically setting up the HMM according to the requirements for the R package [depmixS4](https://cran.r-project.org/web/packages/depmixS4/index.html). See <code>fit_models.R</code> for examples of how to run <code>global_hmm.R</code> with both provided initial parameters and randomized initial parameters.
