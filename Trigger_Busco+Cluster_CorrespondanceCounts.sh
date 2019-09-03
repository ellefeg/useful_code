#! /bin/bash

# script: RunClusterBuscoValidation.sh
# usage: nohup RunClusterBuscoValidation.sh {sampleinput} &
# where {sampleinput} is a tab-delimited file of | sampleID | busco.fnodes dir | clustering.fnodes dir |

sample_input=$1

while read sample busco clustering
do
/home/laura/scripts/Busco+Cluster_CorrespondanceCounts.sh "$sample" "$busco" "$clustering"
done < "$sample_input"
