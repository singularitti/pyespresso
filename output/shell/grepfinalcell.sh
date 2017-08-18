#!/usr/bin/env bash 
# Helpful guide: https://github.com/koalaman/shellcheck/wiki/SC1065

folders=()

for folder in vc_*; do
	[[ -d $folder ]] || break # handle the case of no vc_* folders
	folders+=($folder)
done

for folder in "${folders[@]}"; do
	cd "$folder" && echo "$folder"
    awk "/Begin final coordinates/,/End final coordinates/" *.out
    cd ..
done

function red() {
    printf ("%s%s%s", "\033[1;31m", "$1", "\033[0m ")
}

function green() {
    printf ("%s%s%s","\033[1;32m", "$1", "\033[0m ")
}

function blue() {
    printf ("%s%s%s","\033[1;34m", "$1", "\033[0m ")
}