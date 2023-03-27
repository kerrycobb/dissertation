#!/usr/bin/env fish

set name $argv[1]

scp -r kac0070@easley.auburn.edu:anaxyrus/iqtree/$name/\*.contree .
