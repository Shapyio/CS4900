#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=proj02
#SBATCH --output=proj02.txt
#SBATCH --time=30:00

function decrypter {
	y=$(echo "$1 * 4" | bc) # "$1" is to get core # (1-16)
        x=$(echo "$y-3" | bc)
        for i in $(seq "$x" "$y") # Each core # covers different subset of 1-64 in blocks of 4
        do
                for j in {1..128} # Not elegant; but brute forcing last two inputs
                do
                        for k in {1..128}
                        do
                                ./decrypt  $i $j $k w1140sxs
                        done
                done
        done
}

date
echo "Starting decryption..."
for i in {1..16} # For the 16 cores
do 
	decrypter $i &
done
wait
echo "All done"
date
