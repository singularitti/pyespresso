#!/usr/bin/env bash 

folders=()

for folder in k=*; do
	[[ -d $folder ]] || break # handle the case of no k=_* folders
	folders+=($folder)
done

if [ -f total_energy ]; then
	rm total_energy
else
	touch total_energy
fi

for folder in "${folders[@]}"; do
	cd "$folder" || return
	grep -E "!.*total energy\s+=\s+(-?\d+\.\d+)" ./*.out >>../total_energy
	cd ..
done
