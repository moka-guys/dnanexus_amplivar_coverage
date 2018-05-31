#!/usr/bin/python2

'''
Author:         Amy Slater
Created:        April 2017
Modified:       May 2018
Description:    Parses Amplivar output files and VarScan Or Vardict vcf files to consolidate coverage information into a text report.
                Uses a look up file to annotate amplicon with cDNA region and condons covered and local naming"
Usage:			python path/to/Coverage_rpt.py -c path/to/coverage/files/ -v path/to/vcf/files/ -r path/to/amplicon_lookup_file.txt
'''

# Input files required per sample
# * _flanked.txt
# *.varscan.vcf

import os
import fnmatch
import sys
import getopt
import vcf #pyvcf
from pandas import DataFrame


def readdep(df_flank, rpt, flankfilename):
    '''
    Identifies and reports SWIFT amplicons with sub optimal read depths.
    This is written to the coverage report.
    '''
    # write a header to the report 
    rpt.write("AMPLICON COVERAGE\nFile: " + flankfilename + "\n")
    # dataframe for amplicons with <1000 total read depth or < 500 forward or reverse reads. Specify the columns
    df_1000 = DataFrame(columns = ['Gene', 'cDNA', 'Codons','Amplicon_Genomic_Coordinates', 'TotalReads', 'Forward', 'Reverse'])
    # dataframe for amplicons with <2000 total read depth. Specify the columns
    df_2000 = DataFrame(columns = ['Gene', 'cDNA', 'Codons','Amplicon_Genomic_Coordinates', 'TotalReads', 'Forward', 'Reverse'])
    
    # for each amplicon in the flanking file dataframe
    for index, rd in df_flank.iterrows():
        
        # rd['gene'] is in the format 001_GENE_cDNA_1-100_NM_1234_1-33_chr1:1000-1100
        # split this string to extract the required information
        gene_symbol = rd['gene'].split("_")[1]
        cDNA = rd['gene'].split("_")[2] + "_" + rd['gene'].split("_")[3]
        codons = rd['gene'].split("_")[4] + "_" + rd['gene'].split("_")[5] + "_" + rd['gene'].split("_")[6]
        genomic_coord = rd['gene'].split("_")[7]
        
        # identify amplicons with <1000 total read depth or < 500 forward or reverse reads (rd) = poor covereage
        if rd['total'] < 1000 or rd['forward'] < 500 or rd['reverse'] < 500:    
            # write this to the dataframe
            df_1000.loc[len(df_1000)] = gene_symbol, cDNA, codons, genomic_coord, rd['total'], rd['forward'], rd['reverse']

        # identify amplicon with read depth of < 2000
        elif rd['total'] < 2000:
            # write this to the dataframe
            df_2000.loc[len(df_2000)] = gene_symbol, cDNA, codons, genomic_coord, rd['total'], rd['forward'], rd['reverse']
    
    # convert float columns to integers
    df_1000['TotalReads'] = df_1000['TotalReads'].astype(int)
    df_1000['Forward'] = df_1000['Forward'].astype(int)
    df_1000['Reverse'] = df_1000['Reverse'].astype(int)
    df_2000['TotalReads'] = df_2000['TotalReads'].astype(int)
    df_2000['Forward'] = df_2000['Forward'].astype(int)
    df_2000['Reverse'] = df_2000['Reverse'].astype(int)

    # Write these poorly covered amplicons to file with section headers
    # if there are no poorly covered regions
    if len(df_1000) == 0 and len(df_2000) == 0:
        rpt.write("\nAll amplicons were covered > 2000x total, with > 500x Forward and Reverse\n\n")
    # if there were poorly covered regions
    else:
        # report any from df_1000
        if len(df_1000) > 0:
            # write header
            rpt.write("\nAmplicons covered < 1000x total or with either forward or reverse covered < 500x\n")
            
            # write from dataframe (tab seperated without the index)
            rpt.write(df_1000.to_csv(sep = '\t', index = False))
        # if there aren't any amplicons with coverage between 1000 and 2000 
        if len(df_2000) == 0:
            rpt.write("\nAll other amplicons were covered > 2000x total, with > 500x Forward and Reverse\n\n")
        # report any from df_2000
        elif len(df_2000) > 0:
            # write header
            rpt.write("\nAmplicons with total coverage between 1000 and 2000x\n")
            # write from dataframe (tab seperated without the index)
            rpt.write(df_2000.to_csv(sep = '\t', index = False))
            # state anything not mentioned is covered at at least this level
            rpt.write("\nAll other amplicons were covered > 2000x total, with > 500x Forward and Reverse\n\n")
    


def ampSB(df_flank, rpt, flankfile):
    '''
    Identifies SWIFT amplicons with a F/R strand bias (+/-20%)
    The strand bias is added to the input flanking file.
    '''
    # create dataframe to record amplicons with stand bias (sb). state the columns
    df_sb = DataFrame(columns = ['Gene', 'cDNA', 'Codons','Amplicon_Genomic_Coordinates', 'Forward', 'Reverse', 'StrandBiasRatio'])
    # create dataframe to record amplicons with 0 reads in either F or R 
    df_sb_F_or_R_fail = DataFrame(columns = ['Gene', 'cDNA', 'Codons','Amplicon_Genomic_Coordinates', 'Forward', 'Reverse'])
    # for each row in the flanking dictionary
    for index, flankrow in df_flank.iterrows():
        
        # flankrow['gene'] is in the format 001_GENE_cDNA_1-100_NM_1234_1-33_chr1:1000-1100
        # split this string to extract the required information
        gene_symbol = flankrow['gene'].split("_")[1]
        cDNA = flankrow['gene'].split("_")[2] + "_" + flankrow['gene'].split("_")[3]
        codons = flankrow['gene'].split("_")[4] + "_" + flankrow['gene'].split("_")[5] + "_" + flankrow['gene'].split("_")[6]
        genomic_coord = flankrow['gene'].split("_")[7]
        
        try:
            # identify amplicons with F/R read depth strand bias (sb) (+/-20%) to 4dp
            sb = round(float(flankrow['forward']) / float(flankrow['reverse']), 4)
            # if the strandbias is outside 0.8-1.2
            if sb > 1.2 or sb < 0.8:  # for amplicons with strand bias record cDNA and codons affected
                # Add this to the dataframe
                df_sb.loc[len(df_sb)] = [gene_symbol, cDNA, codons, genomic_coord, flankrow['forward'], flankrow['reverse'], sb]
        
        # if forward or reverse == 0 there will be a zerodivision error
        except ZeroDivisionError:
            # ignore it if both F and R are 0
            if int(flankrow['total']) > 0:
                # otherwise if F or R == 0 add to new df
                df_sb_F_or_R_fail.loc[len(df_sb_F_or_R_fail)] = [gene_symbol, cDNA, codons, genomic_coord, flankrow['forward'], flankrow['reverse']]
        
    # convert counts from floats to integers
    df_sb['Forward'] = df_sb['Forward'].astype(int)
    df_sb['Reverse'] = df_sb['Reverse'].astype(int)
    df_sb_F_or_R_fail['Forward'] = df_sb_F_or_R_fail['Forward'].astype(int)
    df_sb_F_or_R_fail['Reverse'] = df_sb_F_or_R_fail['Reverse'].astype(int)

    # write section header to the report
    rpt.write("\nAMPLICON STRAND BIAS\n0 Forward OR Reverse reads\n")
    # report amplicons failing strand bias
    if len(df_sb_F_or_R_fail) > 0:
        rpt.write("Amplicons where either Forward or Reverse read counts = 0:\n")
        rpt.write(df_sb_F_or_R_fail.to_csv(sep = '\t', index = False))
    else:
        rpt.write("No amplicons had either 0 Forward or Reverse reads\n")

    
    # write section header to the report
    rpt.write("\nAMPLICON STRAND BIAS\nStrand bias defined as +/-20% variation between forward or reverse read count\n")
    # report amplicons failing strand bias
    if len(df_sb) > 0:
        rpt.write("Amplicons with a strand bias +/-20% in forward or reverse\n")
        rpt.write(df_sb.to_csv(sep = '\t', index = False))
        rpt.write("No strand bias observed in other amplicons (+/-20% F/R)\n\n")
    else:
        rpt.write("No strand bias observed in amplicons (+/-20% F/R)\n\n")


def vcf_parse_vcf(vcffile, rpt, flankfile):
    """
    This function reads the VCF file in order to provide allele frequency and strand bias for each variant
    """
       
    # create dataframe for variants failing strand bias check, naming columns
    df_vSB = DataFrame(columns = ['CHROM', 'POS', 'REF', 'ALT', 'StrandBias(FwdRead%)', 'Amplicon'])
    # create dataframe for allele frequency of variants
    df_var = DataFrame(columns = ['Gene', 'cDNA', 'Codons','Amplicon_Genomic_Coordinates', 'CHROM', 'POS', 'REF', 'ALT', 'Alt_Allele_Freq'])
    
    
    # open the flanking file
    with open(flankfile,'r') as flanking:
        # capture the file into a list
        flanking_list=flanking.readlines()
    
    # create dictionary to populate with amplicon names
    amplicon_dict={}
    # parse the list of amplicons 
    for amplicon in flanking_list:
        # skip header
        if amplicon.startswith("gene"):
            pass
        else:
            # capture the amplicon name (first column)
            amplicon_name = amplicon.split("\t")[0]
            # amplicon name is in format 001_GENE_cDNA_1-100_NM_1234_1-33_chr1:1000-11000
            # capture the coordinates as the last element
            amplicon_coord = amplicon_name.split("_")[-1]
            # capture chromosome, start and stop
            chrom = amplicon_coord.split(":")[0]
            start = int(amplicon_coord.split(":")[1].split("-")[0])
            stop = int(amplicon_coord.split(":")[1].split("-")[1])
            # put all this in a dict with the amplicon name as key and dict of chr start and stop
            amplicon_dict[amplicon_name] = {"chrom": chrom, "start" : start, "stop" : stop}
    
    # open the vcf file using vcf package
    myvcf = vcf.Reader(open(vcffile, 'r'))
    
    # for each variant in vcf
    for record in myvcf:
        # pass the record and the lookup table to function which determines the name of overlapping amplicon
        amplicon_name = variant_to_amplicon(record, amplicon_dict)
        # pass record, report, vcffilename,stradbias dataframe and the amplicon name to the vcf_strandbias function which returns the updated dataframe
        df_vSB = vcf_strandbias(record, rpt, vcffile, df_vSB, amplicon_name)
        # pass record, report, vcffilename,allele freq dataframe and the amplicon name to the vcf_freq function which returns the updated dataframe
        df_var = vcf_freq(record, rpt,vcffile, df_var, amplicon_name)
        # convert the POS column in dataframes to integer
        df_vSB['POS'] = df_vSB['POS'].astype(int)
        df_var['POS'] = df_var['POS'].astype(int)
    
    # Write section header
    rpt.write("\nVARIANT STRAND BIAS\nFile: " + vcffile + "\nStrand bias defined as >90% variant supporting reads are from 1 read direction (forward or reverse)\nNOTE: For variants covered by >1 amplicon, each amplicon will be displayed on a separate row.\n")
    # report variants failing strand bias
    if len(df_vSB) > 0:
        # write the variants failing strand bias
        rpt.write(df_vSB.to_csv(sep = '\t', index = False))
        rpt.write("\nNo strand bias observed in other variants.\n")
    else:
        rpt.write("\nNo strand bias observed in variants.\n")

    # report the allele freq of all variants,
    # write header
    rpt.write("\nVARIANT/AMPLICON SUMMARY\nFile:" + vcffile + "\nNOTE: For variants covered by >1 amplicon, each amplicon will be displayed on a separate row.\n")

    # if there any variants...
    if len(df_var) > 0:
        # split variants into on and off target
        df_var_on_target = df_var[df_var['Gene'] != "Off target"]
        df_var_off_target = df_var[df_var['Gene'] == "Off target"]
        
        # print all the on-target variants first
        rpt.write(df_var_on_target.to_csv(sep = '\t', index = False))
        # print all the off-target variants next
        rpt.write(df_var_off_target.to_csv(sep = '\t', index = False))
    
def variant_to_amplicon(record, amplicon_dict):
    """
    This function takes a single variant and identifies any amplicons which overlap the variant
    The name of any overlapping variants are put into a comma seperated string.
    """
    # create empty list to be updated
    Amplivar_Name = []
    # parse the list of flanking  
    for amplicon in amplicon_dict:
        if amplicon_dict[amplicon]["start"] <= int(record.POS) <= amplicon_dict[amplicon]["stop"] and str(record.CHROM) == amplicon_dict[amplicon]["chrom"]:
                Amplivar_Name.append(amplicon)
    
    # return the amplicon name
    return (",").join(Amplivar_Name)


def vcf_strandbias(record, rpt, vcffile, df_vSB,amplicon_name):
    '''
    Recieves a dataframe and a single variant from VCF
    Assesses variant for strand bias (F/R +/- 20%) 
    If strand bias is found adds to dataframe
    '''
    #for each sample in the vcf
    for sample in record.samples:
        # if "varscan" in vcf name 
        if "varscan." in vcffile:
            # find variant allele depth and variant fwd allele depth
            alleledepth_fwd = sample['ADF']
            alleledepth = sample['AD']
        # if "vardict" in vcf name 
        elif "vardict." in vcffile:
            # find variant allele depth and variant fwd allele depth
            alleledepth_fwd = sample['ALD'][0]
            alleledepth = sample['AD'][1]

        # calculate strand bias
        # Calculate the proportion of forward allele reads in the total allele depth
        # convert to a percentage and round to 4dp
        sb = round((float(alleledepth_fwd) / float(alleledepth)) * 100, 4)
        
        # Strand bias is defined as variants with >90% supporting reads in one direction
        if sb >= 90 or sb <= 10:
            # reformat ALT to be a string rather than a list.
            alt = ",".join(str(i) for i in record.ALT)
            # add this to the dataframe df_vSB
            df_vSB.loc[len(df_vSB)] = [record.CHROM, record.POS, record.REF, alt, sb, amplicon_name]
    # return the dataframe
    return df_vSB
        


def vcf_freq(record, rpt, vcffile, df_var, amplicon_name):
    ''' Summarise variants identified, eg which amplicon they are located and the alt allele frequency.'''
    # reformat alt to be a string
    alt = ",".join(str(i) for i in record.ALT)
    # loop through the samples in vcf
    for s in record.samples:
        # If there is an amplicon name for the variant it must be on target
        if amplicon_name:
            # incase variant is covered by > 1 amplicons convert to a list (if only one it'll be a list of 1)
            amplicon_list = amplicon_name.split(",")
            # loop through each amplicon
            for amp_name in amplicon_list:
                # amp_name is in the format 001_GENE_cDNA_1-100_NM_1234_1-33_chr1:1000-1100
                # split this string to extract the required information
                gene_symbol = amp_name.split("_")[1]
                cDNA = amp_name.split("_")[2] + "_" + amp_name.split("_")[3]
                codons = amp_name.split("_")[4] + "_" + amp_name.split("_")[5] + "_" + amp_name.split("_")[6]
                genomic_coord = amp_name.split("_")[7]
                
                # if it's a varscan vcf we need to pull out the FREQ column
                if "varscan." in vcffile:
                    # write to dataframe.
                    df_var.loc[len(df_var)] = [gene_symbol, cDNA, codons, genomic_coord, record.CHROM, record.POS, record.REF, alt, s['FREQ']]
                
                # if it's a vardict vcf we need to pull out the AF column
                elif "vardict." in vcffile:
                    df_var.loc[len(df_var)] = [gene_symbol, cDNA, codons, genomic_coord, record.CHROM, record.POS, record.REF, alt, s['AF']]
        
        # if not on target put blank values 
        else:
            # loop through the samples in vcf
            for s in record.samples:
                # if it's a varscan vcf we need to pull out the FREQ column
                if "varscan." in vcffile:
                    df_var.loc[len(df_var)] = ['Off target', 'N/A', 'N/A', 'N/A', record.CHROM, record.POS, record.REF, alt, s['FREQ']]
                # if it's a vardict vcf we need to pull out the AF column
                elif "vardict." in vcffile:
                    df_var.loc[len(df_var)] = ['Off target', 'N/A', 'N/A', 'N/A', record.CHROM, record.POS, record.REF, alt, s['AF']]
    
    # return the dataframe
    return df_var
    


def main(argv):
    ''' Main function, to locate coverage output files and vcf file for amplivar pipeline'''
    # create empty variables to hold the command line argument values
    covdir = ''
    vcfdir = ''
    
    # check no unrecognised values
    try:
        opts, args = getopt.getopt(argv, "c:v:")
    except getopt.GetoptError:
        print "Unrecognised flag provided"
        sys.exit()
    
    # capture command line values
    for opt, arg in opts:
        if opt == '-c':
            # capture coverage directory
            covdir = arg
        elif opt == '-v':
            # capture vcf directory
            vcfdir = arg

    print covdir
    print vcfdir

    # check vaild directory path has been entered
    assert os.path.isdir(covdir), "Amplivar output directory is not valid"
    # vcf input optional, code will run without it but check (if given) vcfdir is a valid directory
    assert os.path.isdir(vcfdir) or not vcfdir, "vcf input directory is not valid"
    
    # Identify the coverage files provided
    # list of samples
    file_dict = {}
    # loop through all files in coverage directory
    for coverage_file in os.listdir(covdir):
        # capture the sample name by taking everything before _merged_
        sample_name = coverage_file.split("_merged_")[0]
        #coverage file path
        coverage_file_path = os.path.join(covdir,coverage_file)
        # loop through vcf dir to find the associated vcf
        for vcf in os.listdir(vcfdir):
            # if the sample name is in the vcf file name
            if sample_name in vcf:
                # capture full path to vcf
                vcf_file_path = os.path.join(vcfdir,vcf)
        # add the two files for each sample to the dictionary
        file_dict[sample_name] = {"vcf" : vcf_file_path, "flank" : coverage_file_path}

    # report the number of samples
    print "samples identified:"
    print " ,".join(file_dict.keys())

    print 
    # Run coverage report generation for each sample
    for sample_name in file_dict:
        print sample_name
        # create the report file
        rpt = open(str('./' + sample_name + 'coverage_report.txt'), "w")
        print file_dict[sample_name]["flank"]
        print file_dict[sample_name]["vcf"]
        # if the flank file was identified
        if file_dict[sample_name]["flank"]:
            # Dataframe from input flank file - amplicon read depth and amplicon strand bias calling
            df_flank = DataFrame.from_csv(file_dict[sample_name]["flank"], sep = "\t", index_col = None)
            # Call amplicon coverage function
            readdep(df_flank, rpt, file_dict[sample_name]["flank"])
            # Call amplicon strand bias function
            ampSB(df_flank, rpt, file_dict[sample_name]["flank"])
        # if there was no flank file
        else:
            rpt.write("Unable to identify flanked file. Amplicon coverage and strand bias not computed.\n")
            
        # parse VCF file for strand bias and alelle freq
        if file_dict[sample_name]["vcf"]:
            # capture vcf file path
            vcffilepath = file_dict[sample_name]["vcf"]
            # call vcf functions
            vcf_parse_vcf(file_dict[sample_name]["vcf"], rpt, file_dict[sample_name]["flank"])
        # if no vcf
        else:
            rpt.write("Unsuitable number of vcf files found. Variant frequency and strand bias not computed.\n")
        # close the report
        rpt.close()


if __name__ == '__main__':
    main(sys.argv[1:])
