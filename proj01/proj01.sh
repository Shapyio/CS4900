#
# cs4900-01
# Project 01
# Shapiy Sagiev
# w1140sxs
# Due 18 Jan 2024  - accepted up to 10pm Jan19
# System = fry
# Compiler syntax  = NA
# Job Control File = NA
# Additional File  = NA
# Results file     = proj01.txt
#

#!/bin/bash
# VERSION INFO
VERSION=1.0
if [[ $@ == "-v" || $@ == "-version" || $@ == "--v" || $@ == "version" ]]; then
	echo "$0, Version="$VERSION
	exit
fi

# HELP INFO
# Checks if help was called
if [[ $@ == "-h" || $@ == "--h" || $@ == "-help" || $@ == "help" || $@ == "?" || $@ == "-?" ]]; then
	echo "CS 4900"
	echo "Project 01, written by Shapiy Sagiev"
	echo -e "\nTakes in any number of numeric inputs and outputs the sum, average, maximum, and minimum values of the inputted numbers."
	echo -e "\nFlags:"
	echo "  -h --h -help --help ? -?	Opens this useful help menu."
        echo "  -v --v -version --version	Displays program name and version."	
	exit
fi
# Checks if all inputs are numeric
for x in "$@"; do
	if [[ ! $x =~ ^-?[0-9]+$ ]]; then # REGEX was obtained from itslinuxfoss.com, added "-?" for negative #s
		echo "ERROR: Not all inputs are integers."
		exit
	fi
done

# HEADER
echo -e "\n$0 $@"
# SUM
eq=$( echo "$@" | tr ' ' '+')
sum=$( echo "$eq" | bc)
echo "sum="$sum
# AVERAGE
rounded=$(printf "%.1f\n" "$(echo "$sum / $#" | bc -l)") # Rounded average to tenths
if [[ "$rounded" == *".0" ]]; then # If #.0, then output integer
	echo "average="$( echo "$sum / $#" | bc )
else # If #.#, then output decimal
	echo "average="$rounded
fi
# MAX
num=$1
for var in "$@"; do	
	if [[ "$var" -gt "$num" ]]; then
		num=$var
	fi
done
echo "max="$num
# MIN
for var in "$@"; do
	if [ "$var" -lt "$num" ]; then
		num=$var
	fi
done
echo "min="$num
