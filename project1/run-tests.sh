#!/bin/sh
max_threads=10
times_to_execute=10

# Saving script start time
date >> execTime

# Compile program
make

# Remove any previous out directory
rm -r ./out

# Create output directory
mkdir ./out

# Run for each matrix matix
for j in 250 500 1000 1500 2000 3000 4000
do

	for i in $(seq 1 $max_threads)
	do
		echo "Executing for $j Matix - $i Threads"
		in_file_name="./doc/inputs/"$j".txt"
		out_file_name="./out/"$j"-"$i"t.txt"
		echo "./bin/prog $in_file_name $i $times_to_execute $out_file_name"
		./bin/prog $in_file_name $i $times_to_execute $out_file_name
	done
done

date >> execTime

echo " ------------------------------ " >> execTime