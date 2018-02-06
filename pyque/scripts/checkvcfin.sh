#!/usr/bin/env bash
# This script checks whether the vc-relax output gives correct result,
# by examining if there is "End final coordinates" in the output.
# Helpful guides: https://github.com/koalaman/shellcheck/wiki/SC2045,
# http://jafrog.com/2013/11/23/colors-in-terminal.html,
# and https://github.com/koalaman/shellcheck/wiki/SC2068

folders=()

for folder in vc_*; do
	[[ -d $folder ]] || break # handle the case of no vc_* folders
	folders+=($folder)
done

for folder in "${folders[@]}"; do
	cd "$folder" || return
	echo -e "\033[34;1m" "Folder: $folder" "\033[0m"
	if grep -q "End final coordinates" ./*.out; then
		echo -e "\033[32;1m" "Found \"End final coordinates\", finished!" "\033[0m"
	elif grep -q "Maximum number of iterations reached, stopping" ./*.out; then
		echo -e "\033[35;1m" "Maximum number of iterations reached!" "\033[0m"
	elif [ -f CRASH ]; then
		echo -e "\033[38;5;198m" "CRASH file found!" "\033[0m"
	fi
	cd ..
done
