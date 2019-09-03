#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 4 July 2017
# Title: tab2fa.sh
# Goal: To convert a tabular sequence file (from fasta_formatter -t) back to a fasta file
# Usage: ./tab2fa.sh <infile> <out>
# ------------------------------------------------------------------

echo "this script will ONLY work if there is a tab between the sequence header and the sequence"
echo "this should be the case if you have used fasta_formatter -t"
echo "but if you have subsequently used "join" it will not work"
echo "please check your output"

################
###  TAB2FA  ###
################
infile=$1
outfile=$2

date
sed 's/^/>/g' "$infile" | sed 's/\t/\n/g' > "$outfile"
date
