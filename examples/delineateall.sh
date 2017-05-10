#!/bin/bash
# @file delineateall.sh
# @author Jose Vicente <jose.vicenteruiz@fda.hhs.gov>
#
# @version 1.0
#
# @section LICENSE
# This example code is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal Public Domain Dedication. This example is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See DISCLAIMER section below, https://github.com/FDA/ecglib/, https://creativecommons.org/publicdomain/zero/1.0/ and https://www.gnu.org/licenses/gpl-faq.html for more details.
#
# @section DISCLAIMER
# This software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).
# FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
#
# @section DESCRIPTION
# this script parses the annotations list output produced by getdbannotations together with the study clinical data file (SCR-002.Clinical.Data.csv) and calls twavedelineatorphysionet with and without the filtering enabled for each median ECG record.

OLDIFS=$IFS
IFS=","
header="--printheader 1"
while read record egrefid pon qon qoff tpeak tppeak toff
do
    # Assuming EGREFID is incolumn 1 and RR is located in column 21 in the clinical dataset csv file
    IFS=$OLDIFS
    grepcmd="csvtool col 1,21 SCR-002.Clinical.Data.csv | grep "$egrefid""
    rrline=$(eval $grepcmd)
    rr=$(echo "$rrline" | cut -f2 -d,)
    if [ "$rr" != "RR" ]
    then
	# Call T-wave delineator with NO filter
        delineatecmd="./twavedelineatorphysionet --record $record --qon $qon --rpeak 350 --qoff $qoff --rr $rr $header --filterecg 0"
	strout1=$(eval $delineatecmd)
	header="--printheader 0"
	# Call T-wave delineator with filter enabled
	delineatecmd="./twavedelineatorphysionet --record $record --qon $qon --rpeak 350 --qoff $qoff --rr $rr $header --filterecg 1"
	strout2=$(eval $delineatecmd)
        echo -e "$strout1"
        echo -e "$strout2"
    fi
    IFS=","
 done < $1 
 IFS=$OLDIFS
