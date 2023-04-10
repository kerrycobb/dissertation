#!/usr/bin/env python

import fire
import pandas as pd
import os

def make_latex(prefix, group, df):
    df.reset_index(inplace=True)
    df = df.sort_values("Sample ID")
    df["Latitude"] = df["Latitude"].apply(lambda x: '{0:.5f}'.format(x))
    df["Longitude"] = df["Longitude"].apply(lambda x: '{0:.5f}'.format(x))
    last_ix = len(df) - 1
    with open(f"{prefix}-{group}.tex", "w") as fh:
        for ix, row in df.iterrows():
            fh.write("{} & \\textit{{{}}} & {} & {}".format(
                row["Sample ID"], 
                row["Species"], 
                row["Latitude"], 
                row["Longitude"]))
            if ix < last_ix:
                fh.write(" \\\\ \n")

def main(input, output):
    df = pd.read_csv(input)
    df.drop_duplicates(subset="sample_id", inplace=True)
    df["species"] = df["genus"] + " " + df["species"]
    df["prefix"] = df["id"].str.split(" ", n=1, expand=True)[0]
    df = df.rename(columns={
        "id": "Sample ID",
        "species": "Species", 
        "latitude": "Latitude",
        "longitude": "Longitude"})
    kac_df = df[df["prefix"] == "KAC"]
    make_latex(output, "kac", kac_df)
    other_df = df[df["prefix"] != "KAC"]
    make_latex(output, "other", other_df)

if __name__ == '__main__':
    fire.Fire(main)