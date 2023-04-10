#!/usr/bin/env python

import fire
import pandas as pd
import os

def make_latex(path, df):
    df.reset_index(inplace=True)
    df.sort_values("Sample ID", inplace=True)
    df["Latitude"] = df["Latitude"].apply(lambda x: '{0:.5f}'.format(x))
    df["Longitude"] = df["Longitude"].apply(lambda x: '{0:.5f}'.format(x))
    last_ix = len(df) - 1
    with open(path, "w") as fh:
        for ix, row in df.iterrows():
            fh.write("{} & \\textit{{{}}} & {} & {}".format(
                row["Sample ID"], 
                row["Species"], 
                row["Latitude"], 
                row["Longitude"]))
            if ix < last_ix:
                fh.write(" \\\\ \n")

def main(dir):
    dir = dir.strip("/")
    samples = f"{dir}/{dir}.samples.txt"
    kac = f"{dir}/{dir}-kac-samples.tex"
    other = f"{dir}/{dir}-other-samples.tex"
    data = pd.read_csv("~/Desktop/dissertation/toad-data.csv")
    df = pd.read_csv(samples, header=None)
    print(f"Input samples list length: {len(df)}")
    df = df.merge(data, how="left", left_on=0, right_on="sample_id2")
    df.drop_duplicates(subset=0, inplace=True)
    print(f"Output table length: {len(df)}")
    df[["prefix", "num"]] = df["id"].str.split(" ", expand=True)
    df = df.rename(columns={
        "id": "Sample ID",
        "reassign": "Species", 
        "latitude": "Latitude",
        "longitude": "Longitude"})
    include = ["Sample ID", "Species", "Latitude", "Longitude"] 

    kac_df = df[df["prefix"] == "KAC"][include]
    make_latex(kac, kac_df)
    # kac_df.style.to_latex(kac)

    other_df = df[df["prefix"] != "KAC"][include]
    make_latex(other, other_df)
    # other_df.style.to_latex(other)

if __name__ == '__main__':
    fire.Fire(main)