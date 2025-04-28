#!/bin/bash

set -ueo pipefail

OUTPUT_FILE="$1"; shift
touch "$OUTPUT_FILE"

TIME_FORMAT_STRING="%U\t%M"

TEMP_DIR="$(mktemp -d)"
trap '{ rm -rf -- "$TEMP_DIR"; }' EXIT

# Get args

ARGS="$*"
K=""
P=""
COMPLEMENTS=false
while [[ $# -gt 0 ]]; do
    if [[ "$1" == "-k" ]]; then
        shift
        K="$1"
    fi
    if [[ "$1" == "-c" ]]; then
        COMPLEMENTS=true
    fi
    if [[ "$1" == "-p" ]]; then
        shift
        P="$1"
    fi
    shift
done

# Run the right program

if [[ "$ARGS" == *"global"* || "$ARGS" == *"--run-penalty 0"* ]]; then
    # GGMO or LOAC + MO (if rp == 0)
    if [[ "$COMPLEMENTS" = true ]]; then
        COMMAND="./kmercamel $ARGS > "$TEMP_DIR"/ms_raw.txt && ./kmercamel optimize -p "$TEMP_DIR"/ms_raw.txt -a runs -k "$K" -c > "$TEMP_DIR"/ms.txt"
    else
        COMMAND="./kmercamel $ARGS > "$TEMP_DIR"/ms_raw.txt && ./kmercamel optimize -p "$TEMP_DIR"/ms_raw.txt -a runs -k "$K"    > "$TEMP_DIR"/ms.txt"
    fi
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt /bin/sh -c "$COMMAND"
elif [[ "$ARGS" == *"matchtigs"* ]]; then
    # compute matchtigs and turn them into MS
    COMMAND="conda run ggcat build --eulertigs -k "$K" "$P" -j 1 -s 1 -m 64 -p -t "$TEMP_DIR" -o "$TEMP_DIR"/matchtigs.txt && \
            conda run kmercamel spss2ms -k "$K" "$TEMP_DIR"/matchtigs.txt > "$TEMP_DIR"/ms.txt"

    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt /bin/sh -c "$COMMAND"
elif [[ "$ARGS" == *"counting"* ]]; then
    # count distinct kmers
    jellyfish count "$P" -m $K -s 100000 -C -o "$TEMP_DIR"/stats.jf
    SET_SIZE="$(jellyfish stats "$TEMP_DIR"/stats.jf | head -n 2 | tail -n 1 | tr -d '\n' | cut -d: -f2)" # Distinct

    printf "%s\t%s\t:= %d\n" "$(date +"%F:%T")" "$ARGS" "$SET_SIZE" >> "$OUTPUT_FILE"
    exit 0
else
    # LOAC
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/resources.txt \
        ./kmercamel $ARGS > "$TEMP_DIR"/ms.txt
fi

# Compute stats

## Length and runs
LENGTH="$(tail -n +2 "$TEMP_DIR"/ms.txt | wc -m)"
RUNS="$(tail -n +2 "$TEMP_DIR"/ms.txt | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"

## Sizes after compression
tail -n +2 "$TEMP_DIR"/ms.txt > "$TEMP_DIR"/seq.txt

gzip -9 -kq "$TEMP_DIR"/seq.txt; C_GZIP="$(cat "$TEMP_DIR"/seq.txt.gz | wc -c)"
bzip2 -9 -kq "$TEMP_DIR"/seq.txt; C_BZIP2="$(cat "$TEMP_DIR"/seq.txt.bz2 | wc -c)"
xz -9 -kq "$TEMP_DIR"/seq.txt; C_XZ="$(cat "$TEMP_DIR"/seq.txt.xz | wc -c)"
lrzip -Q --zpaq -L9 "$TEMP_DIR"/seq.txt; C_ZPAQ="$(cat "$TEMP_DIR"/seq.txt.lrz | wc -c)"
zstd -22 -q "$TEMP_DIR"/seq.txt; C_ZSTD="$(cat "$TEMP_DIR"/seq.txt.zst | wc -c)"

COMPRESSIBILITY="$(printf "%d %d %d %d %d" "$C_GZIP" "$C_BZIP2" "$C_XZ" "$C_ZPAQ" "$C_ZSTD")"

# Print stats

printf "%s\t%s\t:= %d\t%d\t%s\t%s\t%s\n" "$(date +"%F:%T")" "$ARGS" "$LENGTH" "$RUNS" "$(cat "$TEMP_DIR"/resources.txt)" "$COMPRESSIBILITY" >> "$OUTPUT_FILE"
# printf "%s\t%s\t:= %d\t%d\t%s\t%s\t%s\n" "$(date +"%F:%T")" "$ARGS" "$LENGTH" "$RUNS" "$(cat "$TEMP_DIR"/resources.txt)" >> "$OUTPUT_FILE"

