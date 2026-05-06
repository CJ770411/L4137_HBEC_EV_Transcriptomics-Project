#!/usr/bin/env bash

timestamp() {
    date '+%Y-%m-%d %H:%M:%S: '
}

printf "$(timestamp)"

DIR='asdas'

if [ ! -d "$DIR" ]; then
    printf "\n%sDirectory exists: %s\n\n" "$(timestamp)" "$DIR"
fi