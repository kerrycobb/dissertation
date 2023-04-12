#!/usr/bin/env fish

set name $argv[1]

mkdir $name

scp -r "kac0070@easley.auburn.edu:toad-phyl/dsuite/"$name"/fbranch.svg" $name 