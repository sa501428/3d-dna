#!/bin/bash

#### Description: Scripts forr deep learning based annotation
#### of mismatches. Consists of 2 parts:
####	1. Java code that identifies off-diagonal enrichment at
####       a low resolution, and outputs these regions as
####       individual files to a temp directory. 
####       This code also saves all of the diagonal signal to 
####       the directory. [done]
####	2. Regions are all passed through the neural net and 
####       misjoins are annotated. These are saved to a bed file.
#### Usage: run-deep-mismatch-detector.sh -l <low_res> -r <high_res>
####                   <path-to-hic-file>
#### Dependencies: Juicer_tools; python3 (w/ tensorflow, etc);
#### Input: Juicebox hic file.
#### Parameters: low_res lower resolution on which to find off 
####             diagonal enrichments
####             high_res higher resolution on which to find 
####             and annotate misjoins
#### Output: "Narrow" bed file highlighting mismatch regions 
#### [deep_mismatch_narrow.bed]. These represent candidate regions 
#### of misassembly and might be subject to further filtering.
#### Written by: Muhammad Saad Shamim - shamim@rice.edu. 
#### Version date 01/06/2020.

USAGE="
*****************************************************
This is a wrapper for Hi-C misassembly detection pipeline, version date: Jan 6, 2020. This fragment uses deep learning to generate a mismatch annotation file that will later be overlaid with scaffold boundaries to excise regions spanning misassemblies.

Usage: ./run-mismatch-detector.sh [-h] [-l low_res] [-r high_res] path_to_hic_file

ARGUMENTS:
path_to_hic_file     	Path to .hic file of the current assembly.

OPTIONS:
-h			Shows this help
-l low_res			Sets resolution for the first-pass search of off-diagonal mismatches (default is 25000 bp)
-r high_res		Sets resolution for the precise mismatch localizaton (r<l, default is 1000 bp)

Unprompted
-b NONE/VC/VC_SQRT/KR	Sets which type of contact matrix balancing to use (default KR)

Uses compute-quartile.awk, precompute-depletion-score.awk [[...]] that should be in the same folder as the wrapper script.

*****************************************************
"

## Set defaults
low_res=25000			# default bin size to do a first-pass search for mismatches
high_res=1000	# default bin size to do a second-pass search for mismatches

## Set unprompted defaults
norm="KR"				# use an unbalanced contact matrix for analysis

## HANDLE OPTIONS

while getopts "h:l:r:b:" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
    l)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -l flag was triggered, searching for off-diagonal enrichments at lower resolution of $OPTARG" >&1
            low_res=$OPTARG
        else
            echo ":( Wrong syntax for bin size. Using the default value 25000" >&2
        fi
    ;;
    r)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -r flag was triggered, neural net will receive map regions for annotation at $OPTARG resolution" >&1
            high_res=$OPTARG
        else
            echo ":( Wrong syntax for mismatch localization resolution. Using the default value 1000" >&2
        fi
    ;;
    b)	if [ $OPTARG == NONE ] || [ $OPTARG == VC ] || [ $OPTARG == VC_SQRT ] || [ $OPTARG == KR ]; then
    	    echo ":) -b flag was triggered. Type of norm chosen for the contact matrix is $OPTARG." >&1
			norm=$OPTARG
    	else
    		echo ":( Unrecognized value for -b flag. Running with default parameters (-b NONE)." >&2
    	fi
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## check parameters for consistency
[[ ${low_res} -le ${high_res} ]] && echo >&2 ":( Requested mismatch localization resolution ${high_res} and off diagonal search bin size ${low_res} parameters are incompatible ($low_res < ${high_res}). Exiting!" && echo >&2 "$USAGE" && exit 1

## HANDLE ARGUMENTS: TODO check file format
if [ $# -lt 1 ]; then
    echo ":( Required arguments not found. Please double-check your input!!" >&2
    echo "$USAGE" >&2
    exit 1
fi

hic_file=$1

## CHECK DEPENDENCIES
	type java >/dev/null 2>&1 || { echo >&2 ":( Java is not available, please install/add to path Java to run Juicer and Juicebox. Exiting!"; exit 1; }

###TODO add python check

path_to_scripts=`cd "$( dirname $0)" && pwd`
path_to_vis=$(dirname ${path_to_scripts})"/visualize"
path_to_this_folder=$(dirname ${path_to_scripts})"/deep_edit"

temp_dir="temp_dir"
mkdir ${temp_dir}
juicebox=${path_to_vis}/"juicebox_tools.sh"

if [ ! -f ${juicebox} ] ; then
    echo >&2 ":( Relevant Juicer_tools scripts not found. Exiting!" && exit 1
fi

## DUMP MATRIX FOR ANALYSIS, COMPUTE (1-PCT) QUARTILE AND COMPUTE DEPLETION SCORE ?? substitute for homebrewed dump from mnd? ## TODO: if stick with juicebox_tools ask for proper exit code on fail.
echo "...Initial pass to extract diagonal and regions of interest from .hic file"
${juicebox} shuffle -r ${low_res},${high_res} -k ${norm} -m 500 ${hic_file} null ${temp_dir}
[ $? -ne 0 ] && echo >&2 ":( Initial pass is empty! Perhaps something is wrong with the hic file. Exiting!" && exit 1

python ${path_to_this_folder}/run-deep-mismatch-net-tensorflow.py ${temp_dir} ${high_res} ${path_to_this_folder}/DeepMisjoinNet500x500.h5

rm -r ${temp_dir}
