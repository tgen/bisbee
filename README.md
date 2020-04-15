# bisbee
alternative splicing analysis pipeline

## Prepare input (currently only implemented for spladder)

`python utils/prep.py spladder_counts_file event_type outname spladder_version`

example:

`python utils/prep.py merge_graphs_alt_3prime.counts.hdf5 IR bisbee.IR.counts.txt 2`

## Statistical analysis
### differential splicing

`Rscript stats/diff.R  bisbee_counts_file sample_table outname maxW`

example:

`Rscript stats/diff.R bisbee.IR.counts.txt sample_table.txt bisbee.IR.diff.txt 200`

sample table should have the sample names in the first column and the sample group in the second
