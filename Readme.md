# Dnanexus Amplivar-coverage v1.0

## What does this app do?

- This app runs a file parser and generates a single report per sample summarising the output from amplivar and variant calling into one document.
- This app annotates amplicons with the local 'lab' name, the cDNA regions and the condons the amplicon covers (form henceforth shall be termed *amplicon information*).
- This app highlights amplicons with low coverage or bias in read depth between forward and reverse strands. See below for definitions.
- This app highlights all variants where stand bias is observed. See below for definitions.
- This app summarises all variants identified, allele frequency and annotates with the *amplicon information*.
- This app can take multiple coverage and vcf files as input, generating coverage reports on multiple samples.

## What are typical use cases for this app?

This app generates a single file per sample, collating and summarising the read depth and variant information output as several files by the amplivar pipeline, for testing somatic cancer samples using the SWIFT amplicon panels. 
This app highlights amplicons and variants (if a vcf is supplied) with low read depth , strand bias, which may require further investigation.
Coverage report will only be generated for a sample if a Amplicon read depth tally file (`*flanked.txt`) is supplied. 

## What inputs are required for this app to run?

This app requires the following data:

- Amplicon read depth tally file(s), output from Amplivar (e.g `*_merged_seqprep.fastq.cut.fastq.fna_1_num_grp_flanked.txt`)(REQUIRED).
- vcf file(s) (`*.vcf`) output from either Varscan or VarDict (Optonal)
- Tab delimited look up file detailing the *amplicon information* for the pannel analysed (REQUIRED). File must have the file externtion `*.tsv` and must be in the following format:
	`Chr	Start	End	Gene	Key_Hotspots	Orientation	Swift_LocusID	Amplivar_Name	cDNA	Codons`

**Note**: Amplivar_name MUST match the gene column within the coverage `.txt` files output by Amplivar.


## What does this app output?

The app outputs a `.txt` file per sample summarising:

- Amplicons covered < 1000x total or with either forward or reverse covered < 500x.
- Amplicons with total coverage between 1000 and 2000x.
- Amplicons displaying a strand bias in read depth. Defined as +/-20% variation between forward or reverse read count.
- Variants displaying strand bias. Defined as >90% variant supporting reads are from 1 read direction (forward or reverse).
- Summary of all variants within input vcf, annotated with amplicon name, cDNA region and codons covered by the amplicon and the variant allele frequency. 

## How does this app work?
This app functions as a file parser running python script `Coverage_rpt.py`. This script functions by looping through the coverage and vcf files within input directories and compiles the coverage data and variant data into a single document for each sample. Additionally, this app takes a lookup file denoting the *amplicon information*, and annotates the variant and coverage data. This app requires pandas and pyvcf module to function, these modules are packaged within the app as part of miniconda. 