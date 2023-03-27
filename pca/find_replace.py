#!/usr/bin/env python

import fire
import pandas as pd
import subprocess
import fileinput

def find_replace(file, end=" ", data="~/Desktop/toad-project/toad-data.csv", 
        origCol="sample_id2", replCol="rename"):
    df = pd.read_csv(data)
    df.drop_duplicates(subset=["rename"], inplace=True)
    with open(file, 'r') as f: 
        filedata = f.read()
        for ix, row in df.iterrows():
            find = row[origCol]
            repl = row[replCol]
            filedata = filedata.replace(find + end, repl + end)
    with open(file, 'w') as f:
        f.write(filedata)

if __name__ == "__main__":
    fire.Fire(find_replace)
