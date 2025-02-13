#!/bin/bash

INPUT="$1"
K="$2"
TEST_MODE="$3"
# "C" for counting kmers check - slow

DIR="parameters_testing_outputs"
TIME_FORMAT_STRING="time:\t%E\nmemory:\t%M"

mkdir -p "$DIR";

/usr/bin/time -f "$TIME_FORMAT_STRING" -o "$DIR"/csac_time.txt \
    ./kmercamel -p "$INPUT" -k "$K" -a csac > "$DIR"/csac.txt
/usr/bin/time -f "$TIME_FORMAT_STRING" -o "$DIR"/gg_time.txt \
    ./kmercamel -p "$INPUT" -k "$K" > "$DIR"/gg_raw.txt && \
    ./kmercamel optimize -p "$DIR"/gg_raw.txt -a runs -k "$K" > "$DIR"/gg.txt

echo =================

if [ "$TEST_MODE" = "C" ]; then
    O1="$(./count_noncomplement_kmers.py -p "$DIR"/gg.txt -k "$K")"
    O2="$(./count_noncomplement_kmers.py -p "$DIR"/csac.txt -k "$K" -t)"
fi

L1="$(cat "$DIR"/csac.txt | tail -n 1 | wc -m)"
L2="$(cat "$DIR"/gg.txt | tail -n 1 | wc -m)"

R1="$(cat "$DIR"/csac.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R2="$(cat "$DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"


if [ "$TEST_MODE" = "C" ]; then
    echo kmers:
    printf "% 10d csac\n"   "$O1"
    printf "% 10d gg\n"     "$O2"
fi

echo lengths:
printf "% 10d csac\n"   "$L1"
printf "% 10d gg\n"     "$L2"

echo runs:
printf "% 10d csac\n"   "$R1"
printf "% 10d gg\n"     "$R2"

echo resources - csac:
echo "$(cat "$DIR"/csac_time.txt)"
echo resources - gg:
echo "$(cat "$DIR"/gg_time.txt)"
