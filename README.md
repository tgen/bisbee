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

`Rscript stats/outlierScore.R refSamples.IR.bisbeeFit.txt testSamples.IR.bisbeeCounts.txt  test.ref.IR.bisbeeScores.txt`

Also see jetstream workflow: [workflows/outlier.jst](workflows/outlier.jst)

## Protein Sequence Generation

1. Generate protein sequences

`python prot/build.py bisbee_file event_type aapad outname ensemble_release ref_fasta`

 - bisbee_file: output of utils/prep.py
 - event_type: A3 (alt_3prime), A5 (alt_5prime), ES (exon_skip), IR (intron_retention), or MUT (mutually exclusive exons)
 - aapad: number of amino acids to pad around effected sequence in peptide output
 - outname: name of output file
 - ensemble_release: version of ensemble to use
 - ref_fasta: fasta of reference genome
 
 example: 
 
 `python prot/build.py testSamples.IR.bisbeeCounts.csv IR 9 testSamples.bisbeeProt ensemble_release hg19.fa`
 
 2. find unique protein sequences
 
 `python prot/getUnique.py prot_folder prot_prefix`
  
  - prot_folder: folder with results from prot/build.py
  - prot_prefix: outname used with prot/build.py
  
  example:
  `python prot/getUnique.py bisbeeProt testSamples.bisbeeProt
  
  3. select top transcript effect per event
  
  `python prot/getTop.py prot_folder prot_prefix`
  
  - prot_folder: folder with results from prot/build.py
  - prot_prefix: outname used with prot/build.py
  
  example:
  `python prot/getTop.py bisbeeProt testSamples.bisbeeProt
  
Also see jetstream workflow: [workflows/prot.jst](workflows/prot.jst)

## Filtering and annotation
### differential splicing filtering

`python utils/filtDiff.py diff_folder diff_prefix thresh`

 - diff_folder: folder with differential splicing results
 - diff_prefix: prefix of differential splicing file names
 - thresh: threshold for selecting significant events
 
 example:
 
 `python utils/filtDiff.py bisbeeDiff testSamples 8`
 
 ### outlier results filtering

`python utils/filtOut.py out_folder out_prefix thresh sample_count [sample_file] [select_group] [exclude_group] [exclude_count]`

 - out_folder: folder with outlier results
 - out_prefix: prefix of outlier scores file names
 - thresh: threshold for selecting outliers
 - sample_count: number of samples passing outlier threshold to include event in output
 - sample_file: tab delimited text file with sample ids in first column and sample group in second column, no header (optional)
 - select_group: group of samples to look for outliers, required if sample_file is provided
 - exclude_group: exclude events that have more than "exclude_count" outliers in this group (optional)
 - exclude_count: maximum outliers in exclude group to allow, required if exclude_group is provided
 
 example:
 
 `python utils/filtOut.py bisbeeOut test.ref 8`
 
 #### annotation
 
 `python utils/annotate.py bisbee_out prot_folder prot_prefix`
 
 - bisbee_filt: output file from utils/filtDiff.py or utils/filtOut.py
 - prot_folder: folder with results from prot/getTop.py
 - prot_prefix: outname used with prot/getTop.py
 
 example:
 
 `python utils/annotate.py testSamples.bisbeeDiff.thresh8.csv bisbeeProt testSamples.bisbeeProt`
 
 output file columns:
  - event_cat: type of splicing event
    - Alt: alternate 3 or 5 prime splice site
    - MutEx: mutually exclusive exons
    - ExonInc: exon inclusion in the group indicated in the "group_higher" column
    - ExonSkip: exon skipped in the group indicated in the "group_higher" column
    - IntronRet: intron retained in the group indicated in the "group_higher" column
    - IntronExc: intron excluded in the group indicated in the "group_higher" column
  - group_higher: sample group with more of the isoform resulting in the sequence in the "mutPept" column
  - aa_change_type: type of amino acid change relative to ensembl
    - Canonical: protein coding event with both isoforms are in ensembl
    - Novel: the splice event generates novel transcript and novel amino acid sequence
    - Other: protein coding event that does not fit neatly into novel or canonical criteria
    - ND: effect was not determined
    - None: non-coding transcripts, silent, or protein loss events
  - effect_cat:
    - Substitution: amino acid substitution between isoforms, with the group indicated in the "group_higher" column having more of the isoform with sequence in the "mutPept" column, and the other group having more of the isoform with the sequence in the "wtPept" column
    - Insertion: amino acid insertion/deletion, with the group indicated in the "group_higher" column having more of the isoform with the insertion
    - Deletion: amino acid insertion/deletion, with the group indicated in the "group_higher" column having more of the isoform with the deletion
    - Truncation: amino acid sequence is truncated
    - FrameDisruption: splice event causes use of a different stop codon and altered amino acid sequence
    - ProteinLoss: splice event causes loss of start or stop codon
    - NonCoding: splice event only effects non-coding transcripts
    - NA: effect not available
    
