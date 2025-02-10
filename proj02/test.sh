#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=proj02
#SBATCH --output=test.txt
#SBATCH --time=10:00

function decrypter {
       	start=`date +%s`	
	echo "Core:"$1", Starting at: "$start
	for j in {1..128} # Not elegant; but brute forcing last two inputs
	do
	       	for k in {1..128}
	       	do
	               	./decrypt  $i $j $k w1140sxs
	       	done
	done
	end=`date +%s`
	echo "Core:"$1", Ending at: "$end
	echo "Core:"$1", runtime:"$((end-start))
end
}

echo "Starting decryption..."
for i in {1..16} # For the 16 cores
do
        decrypter $i &
done
wait
echo "All done"
