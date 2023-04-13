#!/usr/bin/env fish

set name $argv[1]

scp -r "kac0070@easley.auburn.edu:toad-phyl/dsuite/"$name"/fbranch.svg" $name-fbranch.svg 
rsvg-convert -f pdf -o $name-fbranch.pdf $name-fbranch.svg