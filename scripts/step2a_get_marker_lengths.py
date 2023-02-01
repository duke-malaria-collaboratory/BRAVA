import os
import pandas as pd

def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

refs = snakemake.params["refs"]
primer_file = open(snakemake.params["primers"][0], "r")
primers = primer_file.readlines()
primers = [s.replace(';optional', '') for s in primers]
primers = [s.split('...') for s in primers]
primArray = [None] * (len(primers) * 2)
j = 0
for i in range(len(primers)):
    if i % 2 != 0:
        for k in primers[i]:
            primArray[j] = k.replace('\n', '')
            rev_comp = rev_compl(primArray[j])
            j += 1
            primArray[j] = rev_comp
            j += 1

dict = {}

if os.path.exists(refs):
    targets = os.listdir(refs)
    for target in targets:
        with open(refs + "/" + target) as f:
            next(f)
            sequence = f.read().replace("\n", "")
            for i in primArray:
                sequence = sequence.replace(i.replace("\n", ""), '')
            dict[target.replace(".fasta", "")] = len(sequence)
else:
    sys.exit("\t**Reads not found**")

df = pd.DataFrame(data=dict, index=[0])
print(df)
df.to_csv(snakemake.output[0], index=False)