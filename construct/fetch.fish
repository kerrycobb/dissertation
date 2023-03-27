#!/usr/bin/env fish

set name $argv[1]

mkdir $name 

scp -r kac0070@easley.auburn.edu:anaxyrus/construct/$name/\*.pdf $name 
scp -r kac0070@easley.auburn.edu:anaxyrus/construct/$name/\*.Robj $name 
scp -r kac0070@easley.auburn.edu:anaxyrus/construct/$name/\*.RData $name 