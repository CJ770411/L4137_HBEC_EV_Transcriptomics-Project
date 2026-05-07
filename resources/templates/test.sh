#!/usr/bin/env bash

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


DIR1='[]'
DIR2='[]'


DIR_LIST=("$DIR1" "$DIR2")

for directory in "${DIR_LIST[@]}"; do
    if [ -d "${directory}" ]; then
        printf "\n%s: Directory exists: %s\n\n" "$(timestamp)" "${directory}"
    else
        mkdir -p "${directory}"
        printf "\n%s - Directory created: %s\n\n" "$(timestamp)" "${directory}"
    fi
done