#!/bin/bash

set -ueo pipefail

OUTPUT_FILE="$1"; shift
touch "$OUTPUT_FILE"

ARGS="$*"

PROGRAM="-"
if [[ "$*" == *"global"* ]]
then
    PROGRAM="G"
fi
# "G" - global greedy

TIME_FORMAT_STRING="%U %M"

TEMP_DIR="$(mktemp -d)"
trap '{ rm -rf -- "$TEMP_DIR"; }' EXIT


if [[ "$PROGRAM" == *"G"* ]]; then
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
        COMMAND="./kmercamel $ARGS > "$TEMP_DIR"/ms_raw.txt && ./kmercamel optimize -p "$TEMP_DIR"/ms_raw.txt -a runs -k "$K" -c > "$TEMP_DIR"/ms.txt"
    else
        COMMAND="./kmercamel $ARGS > "$TEMP_DIR"/ms_raw.txt && ./kmercamel optimize -p "$TEMP_DIR"/ms_raw.txt -a runs -k "$K"    > "$TEMP_DIR"/ms.txt"
    fi
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt /bin/sh -c "$COMMAND"
else
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt \
        ./kmercamel $ARGS > "$TEMP_DIR"/ms.txt
fi

L="$(cat "$TEMP_DIR"/ms.txt | tail -n 1 | wc -m)"
R="$(cat "$TEMP_DIR"/ms.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"

printf "%s\t%s\t:= %d\t%d\t%s\n" "$(date +"%F:%T")" "$ARGS" "$L" "$R" "$(cat "$TEMP_DIR"/resources.txt)" >> "$OUTPUT_FILE"
