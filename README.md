# bisbee
alternative splicing analysis pipeline

## Prepare input (currently only implemented for spladder)

`python utils/prep.py spladder_counts_file event_type outname spladder_version [sample_file] [chunk_num] [chunksize]`

 - spladder_counts_file: \*counts.hdf5 file from spladder
 - event_type: A3 (alt_3prime), A5 (alt_5prime), ES (exon_skip), IR (intron_retention), or MUT (mutually exclusive exons)
 - outname: name for output file
 - spladder_version: 1 for counts files from SplAdder 1.x to correct for changes in later versions
 - sample_file: file with list of samples to include in output if only a subset of the samples from the counts file are desired (optional)
 - chunk_num: output events starting at number chunk_num x chunk_size (optional - intended for analyzing large files in pipelines)
 - chunk_size: number of events to output (optional)

example:

`python utils/prep.py merge_graphs_alt_3prime.counts.hdf5 IR testSamples.IR.bisbeeCounts.txt 2`

## Statistical analysis
### differential splicing

`Rscript stats/diff.R  bisbee_counts_file sample_table outname maxW`
 
 - bisbee_counts_file: output from utils/prep.py
 - sample_table: tab delimited text file with sample ids in first column and sample group in second column, no header
 - outname: output file name
 - maxW: constraint on W parameter, recommended value 200

example:

`Rscript stats/diff.R testSamples.IR.bisbeeCounts.csv sample_table.txt testSamples.IR.bisbeeDiff.txt 200`

Also see jetstream workflow: [workflows/diff.jst](workflows/diff.jst)

### outlier analysis
1. Fit model to reference samples

`Rscript stats/outlierFit.R reference_bisbee_counts_file maxBeta outname`

- reference_bisbee_counts_file: output of utils/prep.py for samples in reference set
- maxBeta: constraint on Beta parameter, recommended value 80
- outname: name of output file

example:

`Rscript stats/outlierFit.R refSamples.IR.bisbeeCounts.txt 80 refSamples.IR.bisbeeFit.txt`

2. Score test samples

`Rscript stats/outlierScore.R bisbee_fit_out test_bisbee_counts outname`

 - bisbee_fit_out: output from stats/outlierFit.R on reference samples
 - test_bisbee_counts: output from utils/prep.py on test samples
 - outname: name of output file

example:

`Rscript stats/outlierScore.R refSamples.IR.bisbeeFit.txt testSamples.IR.bisbeeCounts.txt  test.ref.IRbisbeeScores.txt`

Also see jetstream workflow: [workflows/outlier.jst](workflows/outlier.jst)
