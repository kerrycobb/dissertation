#!/usr/bin/env python

import fire
import pandas as pd
import subprocess
import fileinput

def find_replace(file, end=" ", data="/Users/kerry/Desktop/toad-project/toad-data.csv", 
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


# #!/usr/bin/env python

# import fire
# import pandas as pd
# import subprocess

# def find_replace(file, data="/home/kac0070/anaxyrus/data/sample-data.csv", origCol="sample_id2", replCol="rename"):
#     df = pd.read_csv(data)
#     sedCmds = []
#     for ix, row in df.iterrows():
#         find = row[origCol]
#         repl = row[replCol]
#         sedCmds.append("s/{}/{}/g".format(find, repl))
#     subprocess.check_call("sed -i '{}' {}".format(';'.join(sedCmds), file), shell=True)

# if __name__ == "__main__":
#     fire.Fire(find_replace)