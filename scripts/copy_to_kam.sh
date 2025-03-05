#!/bin/bash

set -ueo pipefail

TARGET="kam"
TARGET_DIR="/home/janci/bakalarka/"
DEFAULT_PARAMS="src scripts Makefile create-version.sh"
# does not include "data" which can be unpractically big and does not change often

ALL_PARAMS="$*"
TO_COPY="${ALL_PARAMS:-"$DEFAULT_PARAMS"}"

COMMAND="cd $TARGET_DIR; nohup make reall"

scp -r $TO_COPY "$TARGET":"$TARGET_DIR"

if [ "$#" -eq 0 ]; then
    ssh "$TARGET" "$COMMAND"
fi
