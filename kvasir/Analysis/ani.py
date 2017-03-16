import os
import pandas as pd
from settings import MONGODB as db
import re
from subprocess import Popen, PIPE

species = db["genes"].distinct("species")
name_dict = {sub(r"[. -\(\)]", "_", x):x for x in species}

l = len(species)
data = {s:[0 for x in range(l)] for s in species}
dm = pd.DataFrame(data, index=species)
dm = pd.read_csv("/Users/KBLaptop/Desktop/ani_matrix.csv", index_col=0)

d = "/Users/KBLaptop/Desktop/ani"
for f1 in os.listdir(d):
    for f2 in os.listdir(d):
        if f1.endswith(".fna") and f2.endswith(".fna"):
            s1 = name_dict[f1[:-4]]
            s2 = name_dict[f2[:-4]]
            if s1 == "Glutamicibacter arilaitensis Re117" or s2 == "Glutamicibacter arilaitensis Re117":
                if s1.split()[0] == s2.split()[0]:
                    print(s1)
                    print(s2)
                    if s1 == s2:
                        dm.loc[s1, s2] = 100.0
                        pass
                    else:
                        dist = Popen(
                            ['ruby', '/Users/KBLaptop/computation/Science/enveomics/Scripts/ani.rb',
                            '--seq1', os.path.join(d, f1),
                            '--seq2', os.path.join(d, f2),
                            '--auto', '--quiet'],
                             stdout=PIPE  # xml output
                        ).communicate()[0]

                        try:
                            print(float(dist.strip()))
                            dm.loc[s1, s2] = float(dist.strip())
                        except ValueError:
                            pass

dm.to_csv("/Users/KBLaptop/Desktop/ani_matrix.csv")
