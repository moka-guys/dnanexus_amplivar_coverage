#!/usr/bin/python2

'''
Author:         Amy Slater
Created:        April 2017
Description:    Parses Amplivar output files to consolidate coverage information into a text report.
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


def readdep(df_flank, rpt, df_bed):
    '''Identifies and reports SWIFT amplicons with read depth coverage <1000 reads or <500 reads F or R'''
    rpt.write("AMPLICON COVERAGE\n")
    # dataframe to record amplicons with poor coverage
    df_1000 = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Codons', 'TotalReads', 'Forward', 'Reverse'])
    df_2000 = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Codons', 'TotalReads', 'Forward', 'Reverse'])
    for index, rd in df_flank.iterrows():
        # identify amplicons with <1000 total read depth and < 500 F and R read depth = poor covereage
        if rd['total'] < 1000 or rd['forward'] < 500 or rd['reverse'] < 500:
            # identify which cDNA and Codons are covered by amplicon from bedfile
            for index, bed in df_bed.iterrows():
                if rd['gene'].strip() == bed['Amplivar_Name']:
                    df_1000.loc[len(df_1000)] = rd['gene'], bed['Gene'], bed['cDNA'], bed['Codons'], rd['total'], rd[
                        'forward'], rd['reverse']

        # identify amplicon with read depth of < 2000
        elif rd['total'] < 2000:
            # identify which cDNA and Codons are covered by amplicon from bedfile
            for index, bed in df_bed.iterrows():
                if rd['gene'].strip() == bed['Amplivar_Name']:
                    df_2000.loc[len(df_2000)] = rd['gene'], bed['Gene'], bed['cDNA'], bed['Codons'], rd['total'], rd[
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


def ampSB(df_flank, rpt, flankfile, df_bed):
    '''Identifies SWIFT amplicons with a F/R strand bias (+/-20%)'''
    rpt.write("\nAMPLICON STRAND BIAS\n")
    rpt.write("(Strand bias defined as +/-20% variation between forward or reverse read count)\n")
    # dataframe to record amplicons with stand bias
    df_sb = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Condons', 'StrandBiasRatio'])
    ratio = []
    for index, sb in df_flank.iterrows():
        # identify amplicons with F/R read depth strand bias (+/-20%)
        s = round(float(sb['forward']) / float(sb['reverse']), 4)
        ratio.append(s)
        if s > 1.2 or s < 0.8: #for amplicons with strand bias record cDNA and codons affected
            for index, bed in df_bed.iterrows():
                if sb['gene'].strip() == bed['Amplivar_Name']:
                    df_sb.loc[len(df_sb)] = [sb['gene'], bed['Gene'], bed['cDNA'], bed['Codons'], s]
    # report amplicons failing strand bias
    if len(df_sb) > 0:
        rpt.write("Amplicons with a strand bias +/-20% in forward or reverse" + "\n")
        # print "Amplicons with a strand bias +/-20% in forward or reverse" + "\n"
        rpt.write(df_sb.to_csv(sep='\t', index=False))
        rpt.write("No strand bias observed in other amplicons (+/-20% F/R)\n\n")
        # print df_sb
    else:
        rpt.write("No strand bias observed in amplicons (+/-20% F/R)\n\n")
    # add strand bias ratio column to raw data in flanked file
    df_flank['strandbiasratio'] = ratio
    df_flank.to_csv(flankfile, sep='\t', index=False)


def varSB(vcffile, rpt, df_bed):
    '''Identifies strand bias (F/R +/- 20%) in variants called'''
    rpt.write("\nVARIANT STRAND BIAS\n")
    rpt.write("(stand bias defined as >90% variant supporting reads are from 1 read direction (forward or reverse))\n")
    # Dataframe for variants failing strand bias check
    df_vSB = DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'StrandBias(FwdRead%)', 'Amplicon'])
    myvcf = vcf.Reader(open(vcffile, 'r'))
    for record in myvcf:
        # record amplicons as a list, so if >1 amplicons cover a variant, all are recorded
        amp = []
        for index, bed in df_bed.iterrows():
            # uses variant position and amplicons start and end location to identify corresponding amplicons
            if record.POS >= bed['Start'] and record.POS <= bed['End']:
                amp.append(bed['Amplivar_Name'])

        for sample in record.samples:
            adf = sample['ADF']
            ad = sample['AD']
            sb = round((float(adf) / float(ad)) * 100, 4)

            # Generate a dataframe of variants displaying strand bias.
            # Strand bias is set as variants with >90% supporting reads in one direction
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


def varFreq(vcffile, rpt,df_bed):
    ''' Summaries variants identified,upon which amplicon they are located and the alt allele frequency.'''
    #print vcffile
    rpt.write("\nVARIANT/AMPLICON SUMMARY\n")
    rpt.write("from:" + vcffile + "\n")
    avcf = vcf.Reader(open(vcffile, 'r'))
    df_var = DataFrame(columns=['Amplicon', 'Gene', 'cDNA', 'Codons','CHROM', 'POS', 'REF', 'ALT', 'Alt_Allele_Freq' ])
    for record in avcf:
        alt = ",".join(str(i) for i in record.ALT) #reformat alt
        # record amplicons as a list, so if >1 amplicons cover a variant, all are recorded
        amp = []
        offtarget = False
        for index, bed in df_bed.iterrows():
            # uses variant position and amplicons start and end location to identify corresponding amplicons
            if record.POS >= bed['Start'] and record.POS <= bed['End']:
                offtarget = True
                for s in record.samples:
                    df_var.loc[len(df_var)] = [bed['Amplivar_Name'],  bed['Gene'], bed['cDNA'],bed['Codons'], record.CHROM, record.POS, record.REF, alt, s['FREQ']]
        if not offtarget:
            for s in record.samples:
                    df_var.loc[len(df_var)] = ['Off target', 'N/A', 'N/A','N/A', record.CHROM, record.POS, record.REF, alt, s['FREQ']]


    rpt.write(df_var.to_csv(sep='\t', index=False))


def main(argv):
    ''' Main function, to locate amplivar output files inorder to report coverage and strand bias'''
    ampoutput = ''
    ampdir = ''
    # check: script has been presented with a directory
    try:
        opts, args = getopt.getopt(argv, "o:a:")
    except getopt.GetoptError:
        print "Unrecognised flag provided"
        sys.exit()
    # sets input directory to that specified as the amplivar output directory
    for opt, arg in opts:
        if opt == '-o':
            ampoutput = arg
        elif opt == '-a':
            ampdir = arg
        else:
            assert "Directory not supplied"  # checks an argument has been supplied
            print "Required files not supplied"

    # check: vaild directory path has been entered
    assert os.path.isdir(ampoutput) == True, "Amplivar output directory is not valid"

    ###################################################################################################################
    dirs = []
    dfiles = {}  # generate a dictionary of vcf and flanked files for each sample
    for root, dirnames, filenames in os.walk(ampoutput):
        dirs.append(dirnames)  # generate lists all directs found. Use to to get list of samples use as dict key

    for i in dirs[0]:  # loop through each sample directory identified to find .varscan.vcf and flanked.txt file
        if ampoutput.endswith('/'):
            spath = (str(ampoutput + i + '/'))
        else:
            spath = (str(ampoutput + '/' + i + '/'))
        dvcf = {}
        dflank = {}
        for root, dirnames, filenames in os.walk(spath):
            c = 0
            for filename in fnmatch.filter(filenames, '*.varscan.vcf'):  # find vcf files
                if len(dvcf) > 0:
                    c += 1
                    dvcf['path' + str(c)] = os.path.join(root, filename)
                else:
                    dvcf['path'] = os.path.join(root, filename)
            for filename in fnmatch.filter(filenames, '*flanked.txt'):  # find flanked files
                if len(dflank) > 0:
                    c += 1
                    dflank['path' + str(c)] = os.path.join(root, filename)
                else:
                    dflank['path'] = os.path.join(root, filename)
        dfiles[i] = {'vcf': dvcf, 'flank': dflank}  # add to dictionary as sample(dir): vcf, flanked

    ###################################################################################################################
    # CHECK: black list all dir with no vcf, and flanked and remove from dictionary
    blklst = []
    for k, v in dfiles.iteritems():
        if len(v['vcf']) == 0 and len(v['flank']) == 0:
            print k, "has no flanked or vcf files attributed. No further analysis will be conducted on sample."
            blklst.append(k)
    # remove key, values from dict that are black listed
    for b in blklst:
        del dfiles[b]
    ###################################################################################################################

    # Load Dataframe of amplicon bedfile req for mapping variants and codons to amplicons.
    bedfile = ampdir + "/bin/universal/SwiftPanel_named.tsv"
    df_bed = DataFrame.from_csv(bedfile, sep="\t", index_col=None)

    # Run coverage report generation for each sample
    for k, v in dfiles.iteritems():
        # generate empty coverage report
        if ampoutput.endswith('/'):
            rpt = open(str(ampoutput + k + '/' + k + '_coverage_report.txt'), "w")
        else:
            rpt = open(str(ampoutput + '/' + k + '/' + k + '_coverage_report.txt'), "w")

        # AMPLICON COVERAGE analysis and checks
        if len(v['flank']) == 1:  # check: 1 flanked file found for sample
            for s in v['flank']:
                if "flanked/" + k in v['flank'][s]:  # check sample ID matches flanked file
                    flankfile = v['flank'][s]
                    try:
                        # Dataframe from input Flanked file - amplicon read depth and amplicon strand bias calling
                        df_flank = DataFrame.from_csv(flankfile, sep="\t", index_col=None)
                        # Call amplicon coverage function
                        readdep(df_flank, rpt, df_bed)
                        # Call amplicon strand bias function
                        ampSB(df_flank, rpt, flankfile, df_bed)
                    except:
                        rpt.write("Unable to load Flanked file, amplicon coverage not reported.\n")
                        pass
                else:
                    rpt.write("Flanked.txt file did not match sample ID. Amplicon coverage not reported.\n")
                    pass
        elif len(v['flank']) == 0 or len(v['flank']) > 1:
            # print k, "Number of flanked files found does not equal 1 please check sample output."
            # print "Unsuitable number of flanked files found. Amplicon coverage and strand bias not computed for sample\n"
            rpt.write("Unsuitable number of flanked files found. Amplicon coverage and strand bias not computed.\n")
            pass

        # VARIANT SUMMARY and STRAND BIAS analysis and checks
        if len(v['vcf']) == 1:  # check: 1 vcf file found for sample
            for s in v['vcf']:
                if k + ".blat.varscan.vcf" in v['vcf'][s]:  # check sample ID matches vcf file
                    vcffile = v['vcf'][s]
                    varSB(vcffile, rpt,df_bed)
                    varFreq(vcffile, rpt, df_bed)
                else:
                    rpt.write("vcf file did not match sample ID. Variant frequency and strand bias not reported.\n")
                    pass
        elif len(v['vcf']) == 0 or len(v['vcf']) > 1:
            rpt.write("Unsuitable number of vcf files found. Variant frequency and strand bias not computed.\n")
            pass
        rpt.close()


if __name__ == '__main__':
    main(sys.argv[1:])