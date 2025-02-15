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

/usr/bin/time -f "$TIME_FORMAT_STRING" -o "$DIR"/csac_c_time.txt \
    ./kmercamel -p "$INPUT" -k "$K" -a csac -c > "$DIR"/csac_c.txt
/usr/bin/time -f "$TIME_FORMAT_STRING" -o "$DIR"/gg_c_time.txt \
    ./kmercamel -p "$INPUT" -k "$K" -c > "$DIR"/gg_c_raw.txt && \
    ./kmercamel optimize -p "$DIR"/gg_c_raw.txt -a runs -k "$K" > "$DIR"/gg_c.txt

echo =================

if [ "$TEST_MODE" = "C" ]; then
    O1="$(./count_noncomplement_kmers.py -p "$DIR"/csac.txt -k "$K" -t)"
    O2="$(./count_noncomplement_kmers.py -p "$DIR"/gg.txt -k "$K")"
    O3="$(./count_kmers.py -p "$DIR"/csac_c.txt -k "$K" -t)"
    O4="$(./count_kmers.py -p "$DIR"/gg_c.txt -k "$K")"
fi

L1="$(cat "$DIR"/csac.txt | tail -n 1 | wc -m)"
L2="$(cat "$DIR"/gg.txt | tail -n 1 | wc -m)"
L3="$(cat "$DIR"/csac_c.txt | tail -n 1 | wc -m)"
L4="$(cat "$DIR"/gg_c.txt | tail -n 1 | wc -m)"

R1="$(cat "$DIR"/csac.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R2="$(cat "$DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R3="$(cat "$DIR"/csac_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R4="$(cat "$DIR"/gg_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"


if [ "$TEST_MODE" = "C" ]; then
    echo kmers:
    printf "% 10d csac\t% 10d (c)\n"   "$O1" "$O3"
    printf "% 10d gg\t% 10d (c)\n"     "$O2" "$O4"
fi

echo lengths:
printf "% 10d csac\t% 10d (c)\n"   "$L1" "$L3"
printf "% 10d gg\t% 10d (c)\n"     "$L2" "$L4"

echo runs:
printf "% 10d csac\t% 10d (c)\n"   "$R1" "$R3"
printf "% 10d gg\t% 10d (c)\n"     "$R2" "$R4"

echo resources - csac:
echo "$(cat "$DIR"/csac_time.txt)"
echo resources - gg:
echo "$(cat "$DIR"/gg_time.txt)"

echo "resources - csac (c)":
echo "$(cat "$DIR"/csac_c_time.txt)"
echo "resources - gg (c)":
echo "$(cat "$DIR"/gg_c_time.txt)"
