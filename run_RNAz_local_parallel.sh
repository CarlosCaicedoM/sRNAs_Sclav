#!/bin/bash
rm -rf results_sRNA/RNAz_out
mkdir results_sRNA/RNAz_out

# Define the function to process each alignment file
process_alignment() {
    local aln_file="$1"
    local base_name=$(basename "$aln_file" .aln)
    echo "Processing alignment $base_name"
    rnazSelectSeqs.pl -i 70 -n 4 "results_sRNA/IGRs_alignments/$base_name.aln" | RNAz --both-strands -p 0.9 > "results_sRNA/RNAz_out/${base_name}.out"
}

# Export the function so that it's available for parallel
export -f process_alignment

# Execute in parallel to process all alignment files
find results_sRNA/IGRs_alignments/ -type f -name "*.aln" | parallel process_alignment {}

