## Pegasus : annotation and prediction of oncogenic gene fusions in RNAseq ##

This software is maintained by Francesco Abate and Sakellarios Zairis.
Please cite our [paper](http://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-014-0097-z) if you use this tool in your analysis.


### Requirements: ###
---

- UNIX like operating system
- java
- perl
- python 2.7.x 
    - numpy
    - pandas
    - scikit-learn
- output reports from one of the following fusion detection tools
    - chimerascan 0.4.5
    - defuse 0.4.3
    - bellerophontes 0.4.0
- human genome and annotation files (hg19)
    - hg19.fa
    - hg19.fa.fai
    - Ensembl Homo_sapiens.GRCh37.60.chr.gtf


### Setup: ###
---

Clone this repository and do not alter the directory structure.
Locally train the classifier as follows:

```bash
$ cd learn
$ python train_model.py
```

This will create two serialized data structures in the current directory, in the form of pickle files.
The wrapper for running the pipeline is pegasus.pl and its command line arguments can be seen by executing the file.
Each run of Pegasus will require a configuration file based on the template provided, as well as the creation of a data specification file.
An example of the data specification file can be found in sample_pipeline_input/data_spec.txt.


### Usage: ###
---

Copy the sample configuration file to the directory of your run and modify the required fields:

1. set the path to the cloned Pegasus repository
2. set the paths for the human genome and annotation reference files
3. set the sample_type to be analyzed by Pegasus (this string is matched to a descriptor field in the data_spec.txt input file)

Construct a properly formatted data specification file for the samples to be analyzed in the run.
A sample invocation of Pegasus from the command line would look like this:

```bash
$ pegasus.pl -c config.txt -d data_spec.txt -l log_folder -o output_folder
```

Pegasus is designed to be interrupted at any stage.
When re-running or re-starting the pipeline a second time it is understood that all previously completed steps will be skipped.
To effect a complete re-run from the first step, remove the contents of the log and output folders and then invoke pegasus.pl.


### Output: ###
---

A successfully completed run will produce as the final output a file called pegasus.output.txt.
The report will contain the "Pegasus driver score" as its first column, along with many attributes and annotations of the fusion candidates.
Some of the key fields are listed below:

- DriverScore: Pegasus Driver Score, the value spans from 0 to 1, where 0 is low oncogenicity and 1 high oncogenicity.
- FusionID: unique fusion identification number.
- Sample_Name: sample name as indicated in the data.txt file.
- Program: fusion detection program as indicated in data.txt file.
- Tot/span_reads: total number of reads encompassing the breakpoint as indicated in the fusion detection tool report.
- Split_reads: total number of split reads detecting the junction breakpoint.
- Chr1: chromosome of five prime gene.
- Chr2: chromosome of three prime gene.
- Gene_Start1: genomic starting position of the five prime gene.
- Gene_End1: genomic ending position of the five prime gene.
- Gene_Start2: genomic starting position of the three prime gene.
- Gene_End2: genomic ending position of the three prime gene.
- Strand1: strand of the five prime gene.
- Strand2: strand of the three prime gene.
- Gene_Name1: gene symbol of the five prime gene.
- Gene_Name2: gene symbol of the three prime gene.
- Gene_Breakpoint1: genomic coordinate of the five prime gene breakpoint.
- Gene_Breakpoint2: genomic coordinate of the three prime gene breakpoint.
- Gene_ID1: official id of the five prime gene.
- Gene_ID2: official id of the three prime gene.
- Sample_Type: sample type as indicated in the data.txt file.
- Sample_occurrency_list: list of sample names where the fusion exactly occurs (same breakpoint coordinates).
- Sample_Type_occurency_list: list of sample types where the fusion exactly occurs (same breakpoint coordinates).
- Kinase_info: information about the presence of a kinase in the fusion (3p_KINASE, 5p_KINASE, BOTH_KINASE, NO_KINASE).
- Transcript_ID1: transcript id of the five prime gene.
- Transcript_ID2: transcript id of the three prime gene.
- Reading_Frame_Info: information about the reading frame. Possible values are in-frame and frame-shifted.
- Protein_Start1: starting position of the fused amino acid sequence of the five prime gene.
- Protein_End1: position of the breakpoint in the amino acid sequence of the five prime gene. 
- Protein_Start2: position of the breakpoint in the amino acid sequence of the three prime gene. 
- Protein_End2: ending position of the fused amino acid sequence of the three prime gene.
- Protein_Sequence: Fused amino acid sequence.
- Exon_Gene1: exon number where the breakpoint of the five prime gene falls.
- Exon_Gene2: exon number where the breakpoint of the three prime gene falls.
- Breakpoint_Region1: region where the breakpoint of the five prime gene falls (3p_UTR, 5p_UTR, CDS, Intron).
- Breakpoint_Region2: region where the breakpoint of the three prime gene falls (3p_UTR, 5p_UTR, CDS, Intron).
- Conserved_Domain1: list of retained domains of five prime gene.
- Lost_Domain1: list of disrupted domains of five prime gene.
- Conserved_Domain2: list of retained domains of three prime gene.
- Lost_Domain2: list of disrupted domains of three prime gene.
