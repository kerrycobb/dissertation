#!/usr/bin/env fish

set name $argv[1]

scp -r "kac0070@easley.auburn.edu:/scratch/phyletica/anaxyrus/pyrad/"$name"_outfiles/"$name".ustr" .