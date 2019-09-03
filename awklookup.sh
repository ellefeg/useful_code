#!/usr/bin/awk -f
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: ?
# Title: awklookup.awk
# Goal: To find a col2 value given the col1 value
# Usage: ./awklookup.sh {key value} {lookupfile}
# ------------------------------------------------------------------

# INTRODUCTION

# I wanted a non-grep solution to look up a value from a two-column file
# with the structure {key \t value} such that if I search the key (col1)
# I return the corresponding col2 value. This is the solution I found on
# Stack Overflow (https://stackoverflow.com/questions/4386185/using-awk-for-a-table-lookup)

# Call it like this:

#awklookup.sh keyval lookupfile
#Example:

#$ cat lookupfile
#aaa     111
#bbb     222
#ccc     333
#ddd     444
#zzz     999
#mmm     888
#$ ./lookup.awk ddd lookupfile
#444
#$ ./lookup.awk zzz lookupfile
#999

# NOTES
# using e.g. 999 as keyvalue doesn't work
# if lookupfile has 3 columns, it will only print the column 2 value
# this is an awk file. DO NOT change the first line of this script (it should be: shbang/usr/bin/awk -f)

BEGIN { key = ARGV[1]; ARGV[1]="" }
$1 == key { print $2 }
