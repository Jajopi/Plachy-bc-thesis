#!/bin/bash

INPUT="$1"
K="$2"
DIR="test_correctness_outputs"
TIME_FORMAT_STRING="time:\t%E\nmemory:\t%M"

mkdir -p "$DIR";

/usr/bin/time -f "$TIME_FORMAT_STRING" -o "$DIR"/gg_time.txt \
    ./kmercamel -p "$INPUT" -k "$K" > "$DIR"/gg_raw.txt && \
    ./kmercamel optimize -p "$DIR"/gg_raw.txt -a runs -k "$K" > "$DIR"/gg.txt
/usr/bin/time -f "$TIME_FORMAT_STRING" -o "$DIR"/csac_time.txt \
    ./kmercamel -p "$INPUT" -k "$K" -a csac > "$DIR"/csac.txt

# ./kmercamel -p "$INPUT" -k "$K" -c > "$DIR"/gg_c.txt
# ./kmercamel -p "$INPUT" -k "$K" -a csac -c > "$DIR"/csac_c.txt

echo =================

O1="$(./count_noncomplement_kmers.py -p "$DIR"/gg.txt -k "$K")"
O2="$(./count_noncomplement_kmers.py -p "$DIR"/csac.txt -k "$K" -t)"
# O3="$(./count_kmers.py -p "$DIR"/gg_c.txt -k "$K")"
# O4="$(./count_kmers.py -p "$DIR"/csac_c.txt -k "$K" -t)"

L1="$(cat "$DIR"/gg.txt | tail -n 1 | wc -m)"
L2="$(cat "$DIR"/csac.txt | tail -n 1 | wc -m)"
# L3="$(cat "$DIR"/gg_c.txt | tail -n 1 | wc -m)"
# L4="$(cat "$DIR"/csac_c.txt | tail -n 1 | wc -m)"

R1="$(cat "$DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R2="$(cat "$DIR"/csac.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
# R3="$(cat "$DIR"/gg_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
# R4="$(cat "$DIR"/csac_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"

echo kmers:
echo "$O1" gg
echo "$O2" csac
# echo "$O3" "gg (c)"
# echo "$O4" "csac (c)"
echo lengths:
echo "$L1" gg
echo "$L2" csac
# echo "$L3" "gg (c)"
# echo "$L4" "csac (c)"
echo runs:
echo "$R1" gg
echo "$R2" csac
# echo "$R3" "gg (c)"
# echo "$R4" "csac (c)"
echo resources:
echo "$(cat "$DIR"/gg_time.txt)"
echo "$(cat "$DIR"/csac_time.txt)"
