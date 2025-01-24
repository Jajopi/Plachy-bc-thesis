#!/bin/bash

INPUT="$1"
K="$2"
DIR="test_correctness_outputs"

mkdir -p "$DIR";

./kmercamel -p "$INPUT" -k "$K" > "$DIR"/gg.txt
./kmercamel -p "$INPUT" -k "$K" -a csac > "$DIR"/csac.txt

#./kmercamel -p "$INPUT" -k "$K" -c > "$DIR"/gg_c.txt
#./kmercamel -p "$INPUT" -k "$K" -a csac -c > "$DIR"/csac_c.txt


O1="$(./count_noncomplement_kmers.py -p "$DIR"/gg.txt -k "$K")"
O2="$(./count_noncomplement_kmers.py -p "$DIR"/csac.txt -k "$K" -t)"
#O3="$(./count_kmers.py -p "$DIR"/gg_c.txt -k "$K")"
#O4="$(./count_kmers.py -p "$DIR"/csac_c.txt -k "$K" -t)"

L1="$(cat "$DIR"/gg.txt | tail -n 1 | wc -m)"
L2="$(cat "$DIR"/csac.txt | tail -n 1 | wc -m)"
#L3="$(cat "$DIR"/gg_c.txt | tail -n 1 | wc -m)"
#L4="$(cat "$DIR"/csac_c.txt | tail -n 1 | wc -m)"

R1="$(cat "$DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R2="$(cat "$DIR"/csac.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
#R3="$(cat "$DIR"/gg_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
#R4="$(cat "$DIR"/csac_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"

echo =================
echo kmers:
echo "$O1" gg
echo "$O2" csac
#echo gg "$O3", csac "$O4" "(complements)"
echo lengths:
echo "$L1" gg
echo "$L2" csac
#echo gg "$L3", csac "$L4" "(complements)"
echo runs:
echo "$R1" gg
echo "$R2" csac
#echo gg "$R3", csac "$R4" "(complements)"
