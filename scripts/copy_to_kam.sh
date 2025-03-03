#!/bin/bash

set -ueo pipefail

TARGET="kam:/home/janci/bakalarka/"
DEFAULT_PARAMS="src scripts Makefile create-version.sh"
# does not include "data" which can be unpractically big

ALL_PARAMS="$*"
TO_COPY="${ALL_PARAMS:-"$DEFAULT_PARAMS"}"

scp -r $TO_COPY "$TARGET"
