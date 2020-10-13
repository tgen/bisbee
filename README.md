# Bisbee
alternative splicing analysis package.  

Bisbee perfoms differential splicing analysis, splicing outlier analysis, and splice isoform protein sequence prediction.  Each of these analysis steps may be run independently or may be combined into a pipeline using the workflows below.  Bisbee also provides utility scripts for extracting data, annotating, filtering and summarizing results.  Bisbee has been developed and tested using splice variant detection data from [SplAdder](https://github.com/ratschlab/spladder).

## Overview
![bisbeeWorkflow](https://github.com/tgen/bisbee/blob/master/bisbee.png)

### Citation
Bisbee: A proteomics validated analysis package for detecting differential splicing, identifying splice outliers, and predicting splice event protein effects
Rebecca F. Halperin, Apurva Hegde, Jessica D. Lang, Elizabeth A. Raupach, C4RCD Research Group, Christophe Legendre, Winnie S. Liang, Patricia M. LoRusso, Aleksandar Sekulic, Jeffrey A. Sosman, Jeffrey M. Trent, Sampathkumar Rangasamy, Patrick Pirrotte, Nicholas J. Schork
bioRxiv 2020.08.13.250167; doi: https://doi.org/10.1101/2020.08.13.250167

## Dependencies
 - R (3.5.2 or later)
   - stats4
   - fitdistrplus
   - extraDistr
 - python (3.6.0 or later)
   - pandas
   - numpy
   - pyensembl
   - Bio

## Workflows
Workflows are provided for use in the [jetstream](https://github.com/tgen/jetstream) pipeline framework. The workflows will divide the dataset into chunks and performs each step in parallel on the chunk.  The filter step at the end of the workflows will pull significant events into a single file.  If jetstream is run on a computing cluster with a slurm scheduler, the "--backend slurm" option may be used to submit each task as a slurm job.  If jetstream is run with the (default) local backend, each task will be launched as a process on the local machine.  

To run the full Bisbee workflow:
`jetstream run [--backend slurm] -C config.yaml workflows/full.jst`

The fields below are required in the config.yaml
```yaml
test_samples: [file path string] text file containing table of samples for testing, first column has sample names and second column has group name (diff and outlier)
ref_samples: [file path string] text file containing list of reference samples names for fitting the outlier model (outlier)
chunk_size: [int] number of events to include in each file chunk for parallel processing (all)
maxW: [float] parameter for diff model, recommended value 200 (diff)
maxBeta: [float] parameter for outlier scoring, recommended value 80 (outlier)
diff_thresh: [float] threshold for differential splicing scores, recommended value 8 (diff)
out_thresh: [float] threshold for filtering outlier scores, recommended value 10 (outlier)
out_count: [int] number of samples meeting outlier criteria required (outlier)
testname: [string] name of test samples for use in output file names (outlier and diff)
refname: [string] name of ref samples for use in output file names (outlier)
bisbeePath: [directory path string] path to bisbee installation (all)
ref_spladder_ver: [1 or 2] version of SplAdder used to detect splice variants in reference samples (all)
test_spladder_ver: [1 or 2] version of SplAdder used to detect splice variants in test samples (all)
aapad: [int] number of amino acids to prepend and postpend to altered amino acids in peptide output (prot)
ensemble: [int] version of ensembl to use to get transcript information (prot)
ref: [file path string] reference genome fasta file
group_list: [list of strings] name of groups to compare, for diff splicing will perform all pairwise comparisons by group, for outlier will summarize outlier counts by group (outlier and diff)
event_list: [list] input data by event type
 - name: [A3, A5, ES, IR, or MUT] name of event type (all)
   test_file: [file path string] spladder counts.hdf5 file for the event type above for the test samples (all)
   n_max_test: [int] floor(event count in test_file divided by the chunk_size) (all)
   ref_file: [file path string] spladder counts.hdf5 file for the event type above for the reference samples (outlier)
   n_max_ref: [int] floor(event count in ref_file divided by the chunk_size) (outlier)
```

Workflows for subsets of the analysis are also provided.  The analysis type requiring the parameter are noted in parenthesis above.

## Prepare input (currently only implemented for spladder)

`python utils/prep.py spladder_counts_file event_type outname spladder_version [sample_file] [chunk_num] [chunksize]`

 - spladder_counts_file: \*counts.hdf5 file from spladder
 - event_type: A3 (alt_3prime), A5 (alt_5prime), ES (exon_skip), IR (intron_retention), or MUT (mutually exclusive exons)
 - outname: prefix for name for output file (bisbeeCounts.csv will be added)
 - spladder_version: 1 for counts files from SplAdder 1.x to correct for changes in later versions
 - sample_file: file with list of samples to include in output if only a subset of the samples from the counts file are desired (optional)
 - chunk_num: output events starting at number chunk_num x chunk_size (optional - intended for analyzing large files in pipelines)
 - chunk_size: number of events to output (optional)

example:

`python utils/prep.py merge_graphs_alt_3prime.counts.hdf5 IR testSamples.IR 2`

## Statistical analysis
### differential splicing

`Rscript stats/diff.R  bisbee_counts_file sample_table outname maxW [group1] [group2]`
 
 - bisbee_counts_file: output from utils/prep.py
 - sample_table: tab delimited text file with sample ids in first column and sample group in second column, no header
 - outname: prefix of name for the output file name (bisbeeDiff.csv will be added)
 - maxW: constraint on W parameter, recommended value 200
 - group1: 1st group to use in comparison from sample_table (optional - use if more than two groups in sample_table)
 - group2: 2nd group to use in comparison from sample_table (optional - use if more than two groups in sample_table)

example:

`Rscript stats/diff.R testSamples.IR.bisbeeCounts.csv sample_table.txt testSamples.IR 200`

Also see jetstream workflow: [workflows/diff.jst](workflows/diff.jst)

### outlier analysis
1. Fit model to reference samples

`Rscript stats/outlierFit.R reference_bisbee_counts_file maxBeta outname [sample_file]`

- reference_bisbee_counts_file: output of utils/prep.py for samples in reference set
- maxBeta: constraint on Beta parameter, recommended value 80
- outname: prefix of name for output file (bisbeeFit.csv will be added)
- sample_file: text file with list of samples to use (optional - use if only a subset of samples in the counts file should be used as the reference samples)

example:

`Rscript stats/outlierFit.R refSamples.IR.bisbeeCounts.csv 80 refSamples.IR.bisbeeFit`

2. Score test samples

`Rscript stats/outlierScore.R bisbee_fit_out test_bisbee_counts outname`

 - bisbee_fit_out: output from stats/outlierFit.R on reference samples
 - test_bisbee_counts: output from utils/prep.py on test samples
 - outname: prefix for name of output file (bisbeeOutlier.csv will be added)

example:

`Rscript stats/outlierScore.R refSamples.IR.bisbeeFit.csv testSamples.IR.bisbeeCounts.csv  test.ref.IR`

Also see jetstream workflow: [workflows/outlier.jst](workflows/outlier.jst)

## Protein Sequence Generation

1. Generate protein sequences

`python prot/build.py bisbee_file event_type aapad outname ensemble_release ref_fasta`

 - bisbee_file: output of utils/prep.py or utils/filtOut.py or utils/filtDiff.py
 - event_type: A3 (alt_3prime), A5 (alt_5prime), ES (exon_skip), IR (intron_retention), or MUT (mutually exclusive exons) or ALL (combined event types file)
 - aapad: number of amino acids to pad around effected sequence in peptide output
 - outname: name of output file
 - ensemble_release: version of ensemble to use
 - ref_fasta: fasta of reference genome
 
 example: 
 
 `python prot/build.py testSamples.IR.bisbeeCounts.csv IR 9 testSamples.bisbeeProt 74 hg19.fa`
 
 2. find unique protein sequences
 
 `python prot/getUnique.py prot_folder prot_prefix`
  
  - prot_folder: folder with results from prot/build.py
  - prot_prefix: outname used with prot/build.py
  
  example:
  `python prot/getUnique.py prot testSamples
  
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
 
 `python utils/filtDiff.py diff testSamples 8`
 
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
    - Alt3: alternate 3  prime splice site
    - Alt5: alternate 5 prime splice site
    - MutEx: mutually exclusive exons
    - *group*ExonInc: exon skip/inclusion with *group* having more inclusion
    - *group*IntronRet: intron retained/excluded with *group* having more retention
  - group_higher: sample group with more of the isoform resulting in the sequence in the "altPept" column
  - aa_change_type: type of amino acid change relative to ensembl
    - Canonical: protein coding event with both isoforms are in ensembl
    - Novel: the splice event generates novel transcript and novel amino acid sequence
    - Other: protein coding event that does not fit neatly into novel or canonical criteria
    - ND: effect was not determined
    - None: non-coding transcripts, silent, or protein loss events
  - effect_cat:
    - Substitution: amino acid substitution between isoforms, with the group indicated in the "group_higher" column having more of the isoform with sequence in the "altPept" column, and the other group having more of the isoform with the sequence in the "refPept" column
    - Insertion: amino acid insertion/deletion, with the group indicated in the "group_higher" column having more of the isoform with the insertion
    - Deletion: amino acid insertion/deletion, with the group indicated in the "group_higher" column having more of the isoform with the deletion
    - Truncation: amino acid sequence is truncated
    - FrameDisruption: splice event causes use of a different stop codon and altered amino acid sequence
    - ProteinLoss: splice event causes loss of start or stop codon
    - NonCoding: splice event only effects non-coding transcripts
    - NA: effect not available
    
