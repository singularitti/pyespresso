#!/usr/bin/env bash 

folders=()
for folder in vc_*; do
	[[ -d $folder ]] || break # handle the case of no vc_* folders
	folders+=($folder)
done

for folder in "${folders[@]}"; do
	cd "$folder" && echo "$folder"
    grep "Parallel version (MPI), running on" ./*.out
	grep "PWSCF        :" ./*.out
    cd ..
done
