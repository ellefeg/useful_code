
## Contents

- [sources](#sources)
- [basic awk, grep & sed](#basic-awk-grep--sed)
- [sort, uniq, cut, etc.](#sort-uniq-cut-etc)
- [echo & log files](#echo--log-files)
- [etc](#etc)
- [sources](#sources)
- [awk & sed for bioinformatics](#awk--sed-for-bioinformatics)


## sources
* https://github.com/stephenturner/oneliners/blob/master/README.md

## basic awk, grep & sed

[[back to top](#contents)]

Extract columns 2, 4, and 5 from file

```bash
awk '{print $2,$4,$5}' file
```
    
Rearrange columns 2, 4, 5 in a different order

    awk 'OFS="\t"; print $5, $4, $2}' file


Print any line where 5th column is (or isn't) equal to `abc123`

    # for pattern matches only
    awk '$5 == "abc123"' file
    # for pattern non-matches only
    awk '$5 != "abc123"' file
    
    
Print any line where the 7th column matches (or doesn't) a regular expression

    # for pattern matches only
    awk '$7  ~ /^[a-f]/' file
    # for pattern non-matches only
    awk '$7 !~ /^[a-f]/' file
    
    
Print rows where column 3 is larger than a given value 

    # col 3 > col 5
    awk '$3 > $5' file
    # col 3 > value
    awk '$3 > 100' file
    # col 3 > variable
    awk '#3 > '$i'' file
    
    # size: >, <, >=, <=
    # equals: ==, !=
    # be careful, '$3 == 100' will find columns that EQUAL 100, but '$3 = 100' will CHANGE column 3 to 100!


Sum column 1

    awk '{sum+=$1} END {print sum}' file


Calculate the average of column 2

    awk '{sum+=$1} END {print sum / NR}' file
    # or
    awk '{x+=$2}END{print x/NR}' file


Replace all occurances of `some` with `text` in file

    sed 's/some/text/g' file.txt


Remove leading/trailing whitespace

    # at the start
    sed 's/^[ \t]*//' file # sed -e optional
    # at the end
    sed 's/[ \t]*$//' file
    # at the start and the end
    sed 's/^[ \t]*//;s/[ \t]*$//' file


Delete blank lines in file

    sed '/^$/d' file
    # or
    grep . filename > newfilename


Delete everything after and including a line containing `EndOfUsefulData`:

    sed -n '/EndOfUsefulData/,$!p' file  
    
Print everything except the first line

    awk 'NR>1' file


Print rows 20-80:

    awk 'NR>=20&&NR<=80' input.txt

Calculate the sum of column 2 and 3 and put it at the end of a row:

    awk '{print $0,$2+$3}' file
    
    
Convert a list into a table, where one column contains a subset of data from the other column

    # example AAD3_Trinity_blah --> AAD3 | AAD3_Trinity_blah
    awk '{print $1, $1}' file | sed 's/_.* /\t/g'
    
    
Remove a given column with `cut` (i.e. keep all but column 5)

    cut -f5 --complement file
    
  
Grep a tab

    grep "$(printf '\t')" file


Find start-to-stop amino acid sequences (L|M (start codon) to * (stop codon))

    # make a tab delimited file
    fasta_formatter -t -i somefile.fa -o somefile.fa.tab
    grep "$(printf '\t')L\|$(printf '\t')M" somefile.fa.tab | grep "\*" > somefile.fa.tab_extract
    # convert back to fasta file (must be tab-delimited table)
    ~/scripts/tab2fa.sh somefile.fa.tab_extract somefile.fa_extract

    
Pull out a subset of characters using awk

    # awk '{print substr(wordNmbr, from, to)}' file
    # e.g. given the text "hello world" pull out "hel"
    awk '{print substr($1, 1, 3)}' file

Find files containing text (`-l` outputs only the file names, `-i` ignores the case `-r` descends into subdirectories)

    grep -lir "some text" *
    

Print a specific line (e.g. line 42) from a file:

    sed -n 42p <file>
    

## sort, uniq, cut, etc.

[[back to top](#contents)]

Count the number of unique lines in file

    cat file | sort -u | wc -l


Find lines shared by 2 files (assumes lines within file1 and file2 are unique; pipe to `wd -l` to count the _number_ of lines shared):

    sort file1 file2 | uniq -d
    
    # Safer
    sort -u file1 > a
    sort -u file2 > b
    sort a b | uniq -d

    # Use comm
    comm -12 file1 file2
    
    # Use join, after sorting each file
    join file1 file2
    
    # join, keeping unpairable likes of file 1:
    join -a 1 file1 file2
    
    # join, and fill blanks with NULL
    join -e'NULL' -o auto file1 file2


Pick 10 random lines from a file:

    sort -R file | head -n 10
    # or 
    shuf file.txt | head -n 10


## echo & log files

[[back to top](#contents)]

Write nicely formatted lines to a log (or other) file

    #! /bin/bash
    somecommand
    cat >> mylog.txt <<COMMENT # can be >> or > depending on write/append
    this is a comment or something to skip
    you can use variables like this $v
    but to combine variables and normal text $(echo $v)_dothis
    note that "these are literal quote marks"
    COMMENT
    morecommands
    

Use newlines correctly in echo

    echo -e hello\nworld #this is the correct one
        # hello
        # world

    echo hello\nworld #not this
        # hello\nworld

    echo -e hello"\n"world #not this
        # hello"
        # "world
        
    echo hello"\n"world #not this
        # hello"\n"world


Expand variables correctly in echo (it only works if the variable is in ""s)

    myvar = $(echo -e hello\nworld)
    echo "$myvar" #this is the correct one
        # hello
        # world
    echo $myvar #not this
        # helloworld


## etc

[[back to top](#contents)]

Avoid multiple redirects (ShellCheck SC2129)

    # instead of:
    echo blah >> file
    date >> file
    cat stuff >> file
    # try
    {
      echo foo
      date
      cat stuff
    } >> file


Create a script of the last executed command:

    echo "!!" > foo.sh
    
    
List or delete all files in a folder that don't match a certain file extension (e.g., list things that are _not_ compressed; remove anything that is _not_ a `.foo` or `.bar` file):

    ls !(*.gz)
    rm !(*.foo|*.bar)
    
    
Number each line in file.txt:

    cat -n file.txt

Don't use ls | grep (ShellCheck SC2010)

    # instead of:
    ls /directory | grep mystring
    # use
    echo /directory/*mystring*

## sources
(1) https://github.com/stephenturner/oneliners/blob/master/README.md

## awk & sed for bioinformatics

Keep only top bit scores in blast hits (best bit score only): (Source 1)

    awk '{ if(!x[$1]++) {print $0; bitscore=($14-1)} else { if($14>bitscore) print $0} }' blastout.txt


Keep only top bit scores in blast hits (5 less than the top): (Source 1)

    awk '{ if(!x[$1]++) {print $0; bitscore=($14-6)} else { if($14>bitscore) print $0} }' blastout.txt


Split a multi-FASTA file into individual FASTA files: (Source 1)

    awk '/^>/{s=++d".fa"} {print > s}' multi.fa

Output sequence name and its length for every sequence within a fasta file: (Source 1)

    cat file.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'



Find start-to-stop amino acid sequences (L|M (start codon) to * (stop codon))

    # make a tab delimited file
    fasta_formatter -t -i somefile.fa -o somefile.fa.tab
    grep "$(printf '\t')L\|$(printf '\t')M" somefile.fa.tab | grep "\*" > somefile.fa.tab_extract
    # convert back to fasta file (must be tab-delimited table)
    ~/scripts/tab2fa.sh somefile.fa.tab_extract somefile.fa_extract
