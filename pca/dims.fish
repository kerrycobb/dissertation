#!/usr/bin/env fish

set input $argv[1]

math (wc -l < $input) / 2 
head -n 1 $input | awk '{print NF}'
