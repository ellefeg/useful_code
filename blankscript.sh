#!/bin/bash
# ------------------------------------------------------------------
# Author: 
# Date: 
# Title: 
# Goal: 
# Usage: 
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------
display_usage() {
	echo "ERROR: 1 argument expected for BLAHSCRIPT - exiting!"

	echo -e "\nUsage:\nsomescript.sh somefile.txt\n"
	}

# Check if correct number of arguments (n = 1) provided
if [ $# != 1 ]; then
	display_usage
    exit 1
fi

somevar=$1
