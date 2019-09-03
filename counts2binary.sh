#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 17 April 2018
# Title: counts2binary.sh
# Goal: To take a tab-delim table and convert all numbers > 1 to 1 - so you have a binary presence-absence matrix
# Usage: counts2binary.sh <my.file>
# ------------------------------------------------------------------

# Renaming variables
myfile=$1

sed 's/\t[1-9][0-9]*/\t1/g' $myfile
