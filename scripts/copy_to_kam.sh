#!/bin/bash

set -ueo pipefail

TARGET="kam"
TARGET_DIR="/home/janci/bakalarka/"
DEFAULT_PARAMS="src scripts Makefile create-version.sh compare_inputs.txt compare_run_penalties.txt verify.py convert_superstring.py"
# does not include "data" which can be unpractically big and does not change often

ALL_PARAMS="$*"
TO_COPY="${ALL_PARAMS:-"$DEFAULT_PARAMS"}"
SKIP_BUILD="${1:-""}"
if [[ "$SKIP_BUILD" == "-" ]]; then
    TO_COPY="$DEFAULT_PARAMS"
fi

COMMAND="cd $TARGET_DIR; nohup make reall"

scp -r $TO_COPY "$TARGET":"$TARGET_DIR"

if [ "$#" -eq 0 ]; then
    ssh "$TARGET" "$COMMAND"
fi
