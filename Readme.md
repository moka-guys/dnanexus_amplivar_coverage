# Dnanexus Amplivar-coverage v1.2

## What does this app do?
- This app runs a Python script which parses output from amplivar and variant calling steps into one document per sample. The script:
  - annotates amplicons with the local 'lab' name, the cDNA regions and the codons the amplicon covers (henceforth shall be termed *amplicon information*).
  - highlights amplicons with low coverage or bias in read depth between forward and reverse strands. See below for definitions.
  - highlights all variants where stand bias is observed. See below for definitions.
  - summarises all variants identified, allele frequency and annotates with the *amplicon information*.

## What are typical use cases for this app?
This app is used on somatic cancer samples being assessed using the SWIFT amplicon panels (processed using amplivar and varscan/vardict).
The app collates and summarises read depth and variant information output from several files created by the amplivar pipeline. 
This app highlights amplicons and variants with low read depth or strand bias, which may require further investigation.


## What inputs are required for this app to run?
This app requires the following data:

- Amplicon read depth tally file(s), output from Amplivar (e.g `*_merged_seqprep.fastq.cut.fastq.fna_1_num_grp_flanked.txt`)
- Amplicon must be named in the format: 001_GENE_cDNA_1-100_NM_1234_1-33_chr1:1000-1100 - This naming is created by amplivar using amplicon naming in the flanking file
- vcf file(s) (`*.vcf`) output from either Varscan or VarDict 

## What does this app output?
The app outputs a `.txt` file per sample summarising:

- Amplicons covered < 1000x total or with either forward or reverse covered < 500x.
- Amplicons with total coverage between 1000 and 2000x.
- Amplicons where either the forward or reverse read count is 0
- Amplicons displaying a strand bias in read depth. Defined as +/-20% variation between forward or reverse read count.
- Variants displaying strand bias. Defined as >90% variant supporting reads are from 1 read direction (forward or reverse).
- Summary of all variants within input vcf, annotated with amplicon name, cDNA region and codons covered by the amplicon and the variant allele frequency. 

## How does this app work?
This app functions as a file parser running python script `Coverage_rpt.py`. This script functions by looping through the coverage and vcf files within input directories and compiles the coverage data and variant data into a single document for each sample. Additionally, this app takes a lookup file denoting the *amplicon information*, and annotates the variant and coverage data. This app requires pandas and pyvcf module to function, these modules are packaged within the app as part of miniconda. 
