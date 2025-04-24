#!/bin/bash

set -ueo pipefail

OUTPUT_FILE="$1"; shift
touch "$OUTPUT_FILE"

ARGS="$*"

TIME_FORMAT_STRING="%U\t%M"

TEMP_DIR="$(mktemp -d)"
trap '{ rm -rf -- "$TEMP_DIR"; }' EXIT


if [[ "$*" == *"global"* || "$*" == *"--run-penalty 0"* ]]; then
    # GGMO or LOAC + MO (if rp == 0)
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
    # if [[ false ]]; then
        # echo TODO compute matchtigs and turn them into MS
    # else
        # LOAC
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt \
            ./kmercamel $ARGS > "$TEMP_DIR"/ms.txt
    # fi
fi

LENGTH="$(tail -n 1 "$TEMP_DIR"/ms.txt | wc -m)"
RUNS="$(tail -n 1 "$TEMP_DIR"/ms.txt | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"

# printf "%s\t%s\t:= %d\t%d\t%s\n" "$(date +"%F:%T")" "$ARGS" "$LENGTH" "$RUNS" "$(cat "$TEMP_DIR"/resources.txt)" >> "$OUTPUT_FILE"

tail -n 1 "$TEMP_DIR"/ms.txt > "$TEMP_DIR"/seq.txt

gzip -9 -kq "$TEMP_DIR"/seq.txt; C_GZIP="$(cat "$TEMP_DIR"/seq.txt.gz | wc -m)"
bzip2 -9 -kq "$TEMP_DIR"/seq.txt; C_BZIP2="$(cat "$TEMP_DIR"/seq.txt.bz2 | wc -m)"
xz -9 -kq -z "$TEMP_DIR"/seq.txt; C_XZ="$(cat "$TEMP_DIR"/seq.txt.xz | wc -m)"
lrzip -Q --zpaq -L8 "$TEMP_DIR"/seq.txt; C_ZPAQ="$(cat "$TEMP_DIR"/seq.txt.lrz | wc -m)"

COMPRESSIBILITY="$(printf "%d\t%d\t%d\t%d" "$C_GZIP" "$C_BZIP2" "$C_XZ" "$C_ZPAQ")"

printf "%s\t%s\t:= %d\t%d\t%s\t%s\n" "$(date +"%F:%T")" "$ARGS" "$LENGTH" "$RUNS" "$(cat "$TEMP_DIR"/resources.txt)" "$COMPRESSIBILITY" >> "$OUTPUT_FILE"

