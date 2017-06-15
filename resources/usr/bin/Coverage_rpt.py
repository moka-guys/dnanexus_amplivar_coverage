#!/usr/bin/python2

'''
Author:         Amy Slater
Created:        April 2017
Modified:       June 2017
Description:    Parses Amplivar output files and VarScan Or Vardict vcf files to consolidate coverage information into a text report.
                Uses a look up file to annotate amplicon with cDNA region and condons covered and local naming"
To run:			python path/to/Coverage_rpt.py -c path/to/coverage/files/ -v path/to/vcf/files/ -r path/to/amplicon_lookup_file.txt
'''

# Input files required per sample
# * _flanked.txt
# *.varscan.vcf

import os
import fnmatch
import sys
import getopt
import vcf
from pandas import DataFrame


def readdep(df_flank, rpt, df_lookup, flankfilename):
    '''Identifies and reports SWIFT amplicons with read depth coverage <1000 reads or <500 reads F or R'''
    rpt.write("AMPLICON COVERAGE\n")
    rpt.write("File: " + flankfilename + "\n")
    # dataframe to record amplicons with poor coverage
    df_1000 = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Codons', 'TotalReads', 'Forward', 'Reverse'])
    df_2000 = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Codons', 'TotalReads', 'Forward', 'Reverse'])
    for index, rd in df_flank.iterrows():
        # identify amplicons with <1000 total read depth and < 500 F and R read depth (rd) = poor covereage
        if rd['total'] < 1000 or rd['forward'] < 500 or rd['reverse'] < 500:
            # identify which cDNA and Codons are covered by amplicon from look up table
            for index, lookup in df_lookup.iterrows():
                if rd['gene'].strip() == lookup['Amplivar_Name']:
                    df_1000.loc[len(df_1000)] = rd['gene'], lookup['Gene'], lookup['cDNA'], lookup['Codons'], rd[
                        'total'], rd[
                                                    'forward'], rd['reverse']

        # identify amplicon with read depth of < 2000
        elif rd['total'] < 2000:
            # identify which cDNA and Codons are covered by amplicon from look up table
            for index, lookup in df_lookup.iterrows():
                if rd['gene'].strip() == lookup['Amplivar_Name']:
                    df_2000.loc[len(df_2000)] = rd['gene'], lookup['Gene'], lookup['cDNA'], lookup['Codons'], rd[
                        'total'], rd[
                                                    'forward'], rd['reverse']

    # report amplicons failing RD
    if len(df_1000) > 0:
        rpt.write("Amplicons covered < 1000x total or with either forward or reverse covered < 500x\n")
        rpt.write(df_1000.to_csv(sep='\t', index=False))
        if len(df_2000) == 0:
            rpt.write("All other amplicons were covered > 2000x total, with > 500x Forward and Reverse\n\n")
            # warn of amplicons with lower covereage than expected
    if len(df_2000) > 0:
        rpt.write("\nAmplicons with total coverage between 1000 and 2000x\n")
        rpt.write(df_2000.to_csv(sep='\t', index=False))
        rpt.write("All other amplicons were covered > 2000x total, with > 500x Forward and Reverse\n\n")
    else:
        rpt.write("All amplicons were covered > 1000x total, with > 500x Forward and Reverse\n\n")


def ampSB(df_flank, rpt, flankfile, df_lookup):
    '''Identifies SWIFT amplicons with a F/R strand bias (+/-20%)'''
    rpt.write("\nAMPLICON STRAND BIAS\n")
    rpt.write("(Strand bias defined as +/-20% variation between forward or reverse read count)\n")
    # dataframe to record amplicons with stand bias (sb)
    df_sb = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Condons', 'Forward', 'Reverse', 'StrandBiasRatio'])
    ratio = []
    for index, flankrow in df_flank.iterrows():
        # identify amplicons with F/R read depth strand bias (sb) (+/-20%)
        sb = round(float(flankrow['forward']) / float(flankrow['reverse']), 4)
        ratio.append(sb)
        if sb > 1.2 or sb < 0.8:  # for amplicons with strand bias record cDNA and codons affected
            for index, lookup in df_lookup.iterrows():
                if flankrow['gene'].strip() == lookup['Amplivar_Name']:
                    df_sb.loc[len(df_sb)] = [flankrow['gene'], lookup['Gene'], lookup['cDNA'], lookup['Codons'],
                                             flankrow['forward'], flankrow['reverse'], sb]
    # report amplicons failing strand bias
    if len(df_sb) > 0:
        rpt.write("Amplicons with a strand bias +/-20% in forward or reverse" + "\n")
        rpt.write(df_sb.to_csv(sep='\t', index=False))
        rpt.write("No strand bias observed in other amplicons (+/-20% F/R)\n\n")
        # print df_sb
    else:
        rpt.write("No strand bias observed in amplicons (+/-20% F/R)\n\n")
    # add strand bias ratio column to raw data in flanked file
    df_flank['strandbiasratio'] = ratio
    df_flank.to_csv(flankfile, sep='\t', index=False)


def vcf_strandbias(vcffile, vcffilepath, rpt, df_lookup):
    '''Identifies strand bias (F/R +/- 20%) in variants called'''
    rpt.write("\nVARIANT STRAND BIAS\n")
    rpt.write("File: " + vcffile)
    rpt.write("\nNOTE: For variants covered by >1 amplicon, each amplicon will be displayed on a separate row.\n")
    rpt.write("(stand bias defined as >90% variant supporting reads are from 1 read direction (forward or reverse))\n")
    # Dataframe for variants failing strand bias check
    df_vSB = DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'StrandBias(FwdRead%)', 'Amplicon'])
    myvcf = vcf.Reader(open(vcffilepath, 'r'))

    for record in myvcf:
        # record amplicons as a list, so if >1 amplicons cover a variant, all are recorded
        amp = []
        # uses variant position and amplicons start and end location to name amplicons
        for index, lookup in df_lookup.iterrows():
            if record.POS >= lookup['Start'] and record.POS <= lookup['End']:
                amp.append(lookup['Amplivar_Name'])

        for sample in record.samples:
            # find variant allele depth and variant fwd allele depth from varscan vcf
            if "varscan." in vcffile:
                alleledepth_fwd = sample['ADF']
                alleledepth = sample['AD']
            # find variant allele depth and variant fwd allele depth from vardict vcf
            elif "vardict." in vcffile:
                alleledepth_fwd = sample['ALD'][0]
                alleledepth = sample['AD'][1]
            # Use variant allele depths to calculate strand bias
            sb = round((float(alleledepth_fwd) / float(alleledepth)) * 100, 4)
            # Generate a dataframe of variants displaying strand bias.
            # Strand bias is defined as variants with >90% supporting reads in one direction
            if sb >= 90 or sb <= 10:
                # reformat ALT and amp to be a string rather than a list.
                alt = ",".join(str(i) for i in record.ALT)
                amp = ", ".join(str(i) for i in amp)
                df_vSB.loc[len(df_vSB)] = [record.CHROM, record.POS, record.REF, alt, sb, amp]

    # report variants failing strand bias
    if len(df_vSB) > 0:
        rpt.write(
            "Variants identified with a strand bias. Variant where >90% supporting reads are from 1 read direction (forward or reverse)\n")
        rpt.write(df_vSB.to_csv(sep='\t', index=False))
        rpt.write("No strand bias observed in other variants called.\n")
    else:
        rpt.write("No strand bias observed in variants called.\n")


def vcf_freq(vcffile, vcffilepath, rpt, df_lookup):
    ''' Summaries variants identified,upon which amplicon they are located and the alt allele frequency.'''
    rpt.write("\nVARIANT/AMPLICON SUMMARY\n")
    rpt.write("File:" + vcffile + "\n")
    rpt.write("NOTE: For variants covered by >1 amplicon, each amplicon will be displayed on a separate row.\n")
    avcf = vcf.Reader(open(vcffilepath, 'r'))
    df_var = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Codons', 'CHROM', 'POS', 'REF', 'ALT', 'Alt_Allele_Freq'])
    # loops through vcf variants
    for record in avcf:
        # reformat alt to be a string
        alt = ",".join(str(i) for i in record.ALT)
        ontarget = False
        # uses variants genomic position compared to amplicons start and end location to identify corresponding amplicons
        # loops through lookup table
        for index, lookup in df_lookup.iterrows():
            if record.POS >= lookup['Start'] and record.POS <= lookup['End']:
                ontarget = True
                for s in record.samples:
                    # run loop to pull out info from varscan vcf
                    if "varscan." in vcffile:
                        df_var.loc[len(df_var)] = [lookup['Amplivar_Name'], lookup['Gene'], lookup['cDNA'],
                                                   lookup['Codons'], record.CHROM, record.POS, record.REF, alt,
                                                   s['FREQ']]
                    # run loop to pull out info from VarDict vcf
                    elif "vardict." in vcffile:
                        df_var.loc[len(df_var)] = [lookup['Amplivar_Name'], lookup['Gene'], lookup['cDNA'],
                                                   lookup['Codons'], record.CHROM, record.POS, record.REF, alt, s['AF']]
        if not ontarget:
            for s in record.samples:
                if "varscan." in vcffile:
                    df_var.loc[len(df_var)] = ['Off target', 'N/A', 'N/A', 'N/A', record.CHROM, record.POS, record.REF,
                                               alt, s['FREQ']]
                elif "vardict." in vcffile:
                    df_var.loc[len(df_var)] = ['Off target', 'N/A', 'N/A', 'N/A', record.CHROM, record.POS, record.REF,
                                               alt, s['AF']]
    # write the amplicon annotated variant dataframe to report
    rpt.write(df_var.to_csv(sep='\t', index=False))


def main(argv):
    ''' Main function, to locate coverage output files and vcf file for amplivar pipeline'''
    covdir = ''
    vcfdir = ''
    swiftlookup = ''
    # check: script has been presented with a directory
    try:
        opts, args = getopt.getopt(argv, "c:v:r:")
    except getopt.GetoptError:
        print "Unrecognised flag provided"
        sys.exit()
    # sets input directory to that specified as the amplivar output directory
    for opt, arg in opts:
        if opt == '-c':
            covdir = arg
        elif opt == '-v':
            vcfdir = arg
        elif opt == '-r':
            swiftlookup = arg
        else:
            assert "Directory not supplied"  # checks an argument has been supplied
            print "Required files not supplied"
    print covdir
    print vcfdir
    print swiftlookup
    # check: vaild directory path has been entered
    assert os.path.isdir(covdir) == True, "Amplivar output directory is not valid"
    if len(vcfdir) > 0:  # vcf input optional, code will run without it
        assert os.path.isdir(vcfdir) == True, "vcf input directory is not valid"
    assert os.path.isfile(swiftlookup) == True, "Swift lookup input is not valid"

    ###################################################################################################################
    lcovfilenames = []
    lsample = []
    dsamples = {}  # generate a dictionary of vcf and flanked files for each sample
    for root, dirnames, filenames in os.walk(covdir):
        lcovfilenames.append(filenames)  # generate lists all files output from amplivar. Use to to get list of samples

    for i in lcovfilenames[0]:
        # split file name after sample ID then add the sample id to a new list
        lsample.append(i.split("_merged_")[0])
    print "samples identified:"
    print lsample

    for i in lsample:  # loop through input directories to find coverage file (flanked.txt) and .vcf file
        # initiate a dictionary for the vcf path (dvcf) and a second for the coverage file path (dflank)
        dvcf = {}
        dflank = {}

        # Find vcf files for each sample. Look through vcf input dir
        for root, dirnames, filenames in os.walk(vcfdir):
            vcfname = i + '*.vcf'
            for filename in fnmatch.filter(filenames, vcfname):  # find vcf files for sample
                # add path to vcf dictionary,  as filename:filepath
                dvcf[filename] = os.path.join(root, filename)

        # Find coverage files for each sample. Look through coverage input dir
        for root, dirnames, filenames in os.walk(covdir):
            covname = i + '*flanked.txt'
            for filename in fnmatch.filter(filenames, covname):  # find flanked file
                # add path to coverage dictionary, as filename:filepath
                dflank[filename] = os.path.join(root, filename)

        dsamples[i] = {'vcf': dvcf,
                       'flank': dflank}  # add dictionarys to a big dictionary as sampleID: vcf dictionary, coverage dictionary
        # print dsamples

    ###################################################################################################################
    # CHECK: black list all samples with no vcf, and flanked and remove from dictionary - should not be required in DNAnexus
    blklst = []
    for sample, values in dsamples.iteritems():
        if len(values['vcf']) == 0 and len(values['flank']) == 0:
            print sample, "has no flanked or vcf files attributed. No further analysis will be conducted on sample."
            blklst.append(sample)
    # remove key, values from dict that are black listed
    for b in blklst:
        del dsamples[b]
    ###################################################################################################################

    # Load Dataframe of amplicon lookup table req for mapping variants and codons to amplicons.
    df_lookup = DataFrame.from_csv(swiftlookup, sep="\t", index_col=None)

    # Run coverage report generation for each sample
    for sample, samplefiles in dsamples.iteritems():

        # generate empty coverage report
        rpt = open(str('./' + sample + 'coverage_report.txt'), "w")

        # AMPLICON COVERAGE analysis and checks
        if len(samplefiles['flank']) == 1:  # check: 1 flanked file found for sample
            for sflankfile in samplefiles['flank']:
                # set file name and path as seperate variables inorder to generate dataframs or pass to function.
                flankfilename = sflankfile
                flankfilepath = samplefiles['flank'][sflankfile]
                try:
                    # Dataframe from input Flanked file - amplicon read depth and amplicon strand bias calling
                    df_flank = DataFrame.from_csv(flankfilepath, sep="\t", index_col=None)
                    # Call amplicon coverage function
                    readdep(df_flank, rpt, df_lookup, flankfilename)
                    # Call amplicon strand bias function
                    ampSB(df_flank, rpt, flankfilepath, df_lookup)
                except:
                    rpt.write("Unable to load Flanked.txt coverage file, amplicon coverage not reported.\n")
                    pass
                    # else:
                    # rpt.write("Flanked.txt file did not match sample ID. Amplicon coverage not reported.\n")
                    # pass
        elif len(samplefiles['flank']) == 0 or len(samplefiles['flank']) > 1:
            rpt.write("Unsuitable number of flanked files found. Amplicon coverage and strand bias not computed.\n")
            pass

        # VARIANT SUMMARY and STRAND BIAS analysis and checks
        if len(samplefiles['vcf']) == 1:  # check: 1 vcf file found for sample
            for svcffile in samplefiles['vcf']:
                # set file name and file path as variables
                vcffilepath = samplefiles['vcf'][svcffile]

                # call vcf functions
                vcf_strandbias(svcffile, vcffilepath, rpt, df_lookup)
                vcf_freq(svcffile, vcffilepath, rpt, df_lookup)

        elif len(samplefiles['vcf']) == 0 or len(samplefiles['vcf']) > 1:
            rpt.write("Unsuitable number of vcf files found. Variant frequency and strand bias not computed.\n")
            pass
        rpt.close()


if __name__ == '__main__':
    main(sys.argv[1:])
