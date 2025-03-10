#!/bin/bash

set -ueo pipefail

touch results.txt

ARGS="$*"

PROGRAM="G"
if [[ "$*" == *"csac"* ]]
then
    PROGRAM="C"
fi
# "G" - global greedy
# "C" - csac

OUTPUT_FILE="results.txt"
TIME_FORMAT_STRING="%U %M"

TEMP_DIR="$(mktemp -d)"
trap '{ rm -rf -- "$TEMP_DIR"; }' EXIT

if [[ "$PROGRAM" == *"C"* ]]; then
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt \
        ./kmercamel $ARGS > "$TEMP_DIR"/ms.txt
else
    K=""
    COMPLEMENTS=false
    while [[ $# -gt 0 ]]; do
        if [[ "$1" == *"-k"* ]]; then
            shift
            K="$1"
        fi
        if [[ "$1" == *"-c"* ]]; then
            COMPLEMENTS=true
        fi
        shift
    done

    if [[ "$COMPLEMENTS" = true ]]; then
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt \
            ./kmercamel $ARGS > "$TEMP_DIR"/ms_raw.txt && \
            ./kmercamel optimize -p "$TEMP_DIR"/ms_raw.txt -a runs -k "$K" -c > "$TEMP_DIR"/ms.txt
    else
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt \
            ./kmercamel $ARGS > "$TEMP_DIR"/ms_raw.txt && \
            ./kmercamel optimize -p "$TEMP_DIR"/ms_raw.txt -a runs -k "$K" > "$TEMP_DIR"/ms.txt
    fi
fi

L="$(cat "$TEMP_DIR"/ms.txt | tail -n 1 | wc -m)"
R="$(cat "$TEMP_DIR"/ms.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"

printf "%s %s := %d %d %s\n" "$(date +"%F:%T")" "$ARGS" "$L" "$R" "$(cat "$TEMP_DIR"/resources.txt)" >> "$OUTPUT_FILE"
