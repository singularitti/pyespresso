#!/usr/bin/env bash 
# Helpful guides: https://github.com/koalaman/shellcheck/wiki/SC1065
# and: https://superuser.com/questions/748724/pass-a-large-string-to-grep-instead-of-a-file-name

folders=()
rm finalcell
touch finalcell
echo "# This file is generated in directory $(pwd)" >finalcell

for folder in vc_*; do
	[[ -d $folder ]] || break # handle the case of no vc_* folders
	folders+=($folder)
done

for folder in "${folders[@]}"; do
	cd "$folder" && echo "Current folder is: $folder"
	contents=$(awk "/Begin final coordinates/,/End final coordinates/" ./*.out)
	GREP_COLOR='31;1' grep -E --color=always "new unit-cell volume.*" <<<"$contents"
	grep "new unit-cell volume.*" <<<"$contents" >>../finalcell
	GREP_COLOR='36;1' grep -E --color=always -A 3 "CELL_PARAMETERS.*" <<<"$contents"
	grep -A 3 "CELL_PARAMETERS.*" <<<"$contents" >>../finalcell
	GREP_COLOR='34;1' grep -E --color=always -A 2 "ATOMIC_POSITIONS.*" <<<"$contents"
	grep -A 2 "ATOMIC_POSITIONS.*" <<<"$contents" >>../finalcell
	cd ..
done
