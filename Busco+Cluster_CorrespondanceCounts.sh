#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date:
# Title:
# Goal: This script makes a BUSCO-GeneClustering count matrix showing, for each busco-silix combo:
	# Number of genes with Busco Family X
	# Number of genes in other Busco families
	# Number of genes not in a Busco family
	# Number of other families containing Busco Family X
# And then calculate sensitivity, specificity and jaccard scores for each dataset
# Usage: nohup ./MyScript.sh {sampleID} {busco fnodes} {clustering fnodes} > nohup.out 2>&1&
# NB: BUSCO's `full_table_<sample>.busco.tsv` output file can easily be converted into an *.fnodes style table (BUSCO ID | sequence ID) by cutting the appropriate columns.
# ------------------------------------------------------------------

sampleID=$1
buscofnodes_unsorted=$2
silixfnodes_unsorted=$3

# ------------------------------------------------------------------
# Prepare correspondance file
# ------------------------------------------------------------------

sort -k2,2 "$buscofnodes_unsorted" > $(basename "$buscofnodes_unsorted").sorted
buscofnodes=$(basename "$buscofnodes_unsorted").sorted

sort -k2,2 "$silixfnodes_unsorted" > $(basename "$silixfnodes_unsorted").sorted
silixfnodes=$(basename "$silixfnodes_unsorted").sorted

join -a 1 -1 2 -2 2 "$silixfnodes" "$buscofnodes" | sed 's/ /\t/g' > "$sampleID"_busco+clust.tab

awk '$3 != NULL' "$sampleID"_busco+clust.tab | awk '$2 != NULL' | cut -f 2-3 | sort -u > "$sampleID"_busco+clust.pairs

# ------------------------------------------------------------------
# Make count table
# ------------------------------------------------------------------

touch "$sampleID"_busco+clust.counts1

echo "cluster busco #buscoX #buscoOther #NonBusco #NonClusterBuscoX" >> "$sampleID"_busco+clust.counts1

while read cluster busco
do
grep "$cluster" "$sampleID"_busco+clust.tab > "$sampleID"_"$cluster".busco+clust
echo "$cluster" "$busco" $(grep -c "$busco" "$sampleID"_"$cluster".busco+clust) $(grep "EO" "$sampleID"_"$cluster".busco+clust | grep -c -v "$busco") $(grep -c -v "EO" "$sampleID"_"$cluster".busco+clust) $(grep "$busco" "$sampleID"_busco+clust.tab | grep -c -v "$cluster") >> "$sampleID"_busco+clust.counts1
rm "$sampleID"_"$cluster".busco+clust
done < "$sampleID"_busco+clust.pairs

sed -i 's/ /\t/g' "$sampleID"_busco+clust.counts1

# ------------------------------------------------------------------
# Tidy the results
# ------------------------------------------------------------------

# Get a list of families with 2+ BUSCOs
awk '$4 > 0' "$sampleID"_busco+clust.counts1 > "$sampleID"_2+BuscoClusters

# Keep only the major BUSCO family
sort -nrk3,3 "$sampleID"_busco+clust.counts1 | sort -u -k1,1 > "$sampleID"_busco+clust.counts

# ------------------------------------------------------------------
# Get Statistics
# ------------------------------------------------------------------

# Specificity, sensitivity and Jaccard scores
# Loewenstein, Y., Portugaly, E., Fromer, M., and Linial, M. (2008). Efficient algorithms for accurate hierarchical clustering of huge datasets: tackling the entire protein space. Bioinformatics 24, i41–i49.
# EDIT 29.08.19: See also Salichos, L., and Rokas, A. (2011). Evaluating Ortholog Prediction Algorithms in a Yeast Model Clade. PLoS One 6.
# TruePos = col3; FalsePos = col4+5; FalseNeg=col6

# Specificity = TP / (TP + FP) = $3 / ($3 + $4 + $5)
# Sensitivity = TP / (TP + FN) = $3 / ($3 + $6)
# Jaccard = TP / (TP + FP + FN) = $3 / ($3 + $4 + $5 + $6)

# Get all data
touch "$sampleID"_busco+clust.familyStats
echo "family specificity sensitivity jaccard" >> "$sampleID"_busco+clust.familyStats
grep "EO" "$sampleID"_busco+clust.counts | awk '{print $1, $3 / ($3+$4+$5), $3 / ($3 + $6), $3 / ($3 + $4 + $5 + $6)}' >> "$sampleID"_busco+clust.familyStats
sed -i 's/ /\t/g' "$sampleID"_busco+clust.familyStats


# Get specificity
echo average specificity = $(cut -f 2 "$sampleID"_busco+clust.familyStats | grep -v "specificity" | awk '{sum+=$1} END {print sum / NR}') >> "$sampleID"_busco+clust.AveStats

# Get sensitivity
echo average sensitivity = $(cut -f 3 "$sampleID"_busco+clust.familyStats | grep -v "sensitivity" | awk '{sum+=$1} END {print sum / NR}') >> "$sampleID"_busco+clust.AveStats

# Get jaccard
echo average jaccard = $(cut -f 4 "$sampleID"_busco+clust.familyStats | grep -v "jaccard" | awk '{sum+=$1} END {print sum / NR}') >> "$sampleID"_busco+clust.AveStats

# ------------------------------------------------------------------
# Tidy up
# ------------------------------------------------------------------

rm $(basename "$buscofnodes_unsorted").sorted $(basename "$silixfnodes_unsorted").sorted "$sampleID"_busco+clust.tab "$sampleID"_busco+clust.pairs "$sampleID"_busco+clust.counts1

echo -e "Busco-Clustering analysis complete.\n"$sampleID"_busco+clust.counts = BUSCO-Clustering counts.\n"$sampleID"_busco+clust.familyStats = Row-wise family scores.\n"$sampleID"_busco+clust.AveStats = average Sensitivity, Specificity and Jaccard scores for "$sampleID"."
