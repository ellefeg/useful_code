#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 20 July 2017
# Title: ColumnValues.sh
# Goal: To calculate the average value of a column of data
# Usage: ./ColumnValues.sh <ColNumber> <my.file>
# ------------------------------------------------------------------

# Renaming variables
col=$1
myfile=$2

# Calculate total column sum
echo "Column total:"
cut -f $col $myfile | awk '{sum+=$1} END {print sum}'

# Calculate average column value
echo "Column average:"
cut -f $col $myfile | awk '{sum+=$1} END {print sum / NR}'
