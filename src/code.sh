#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


# Download inputs

dx-download-all-inputs 

#make directory to hold input fastqs 
mkdir coveragein vcfsin

#move all coverage data inputs
for input in /home/dnanexus/in/coverage/*; do if [ -d "$input" ]; then mv $input/* coveragein/; fi; done
#move all vcf inputs if supplied
for input in /home/dnanexus/in/vcfs/*; do if [ -d "$input" ]; then mv $input/* vcfsin/; fi; done
# set -v variable to dir path
if [[ "$vcfs" != "" ]]; then
	allvcfs="-v /home/dnanexus/vcfsin"
else allvcfs=""
fi;

# Add Packaged miniconda to PATH
export PATH="/usr/bin/Miniconda/bin:$PATH"
# test conda
python -m conda list

# run coverage_rpt script 
mark-section "Run file parser"
python /usr/bin/Coverage_rpt.py -c /home/dnanexus/coveragein $allvcfs

## move the output coverage reports and upload
mark-section "Upload output"
#make output folders
mkdir -p out/coverage_rpt/coverage_rpt/ 
# move file
mv ./*coverage_report.txt ~/out/coverage_rpt/coverage_rpt/

dx-upload-all-outputs --parallel
mark-success
