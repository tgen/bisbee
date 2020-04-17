# bisbee
alternative splicing analysis pipeline

## Prepare input (currently only implemented for spladder)

`python utils/prep.py spladder_counts_file event_type outname spladder_version`

example:

`python utils/prep.py merge_graphs_alt_3prime.counts.hdf5 IR testSamples.IR.bisbeeCounts.txt 2`

## Statistical analysis
### differential splicing

`Rscript stats/diff.R  bisbee_counts_file sample_table outname maxW`

example:

`Rscript stats/diff.R testSamples.IR.bisbeeCounts.txt sample_table.txt testSamples.IR.bisbeeDiff.txt 200`

sample table should have the sample names in the first column and the sample group in the second

Also see jetstream workflow: [workflows/diff.jst](workflows/diff.jst)

### outlier analysis
1. Fit model to reference samples

`Rscript stats/outlierFit.R reference_bisbee_counts_file maxBeta outname`

example:

`Rscript stats/outlierFit.R refSamples.IR.bisbeeCounts.txt 80 refSamples.IR.bisbeeFit.txt`

2. Score test samples

`Rscript stats/outlierScore.R refSamples.IR.bisbeeFit.txt testSamples.IR.bisbeeCounts.txt  test.ref.IRbisbeeScores.txt`

Also see jetstream workflow: [workflows/outlier.jst](workflows/outlier.jst)
