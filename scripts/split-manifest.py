from pprint import pprint
from itertools import zip_longest

def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

manifest = 'prostate_cancer_manifest.tsv'
chunk_size= 87 # roughly split into 4 even chunks

with open(manifest) as f:
    lines = [x for x in f
             if 'mini' not in x]
    head = lines[0]
    tail = lines[1:]

    for i, group in enumerate(grouper(chunk_size, tail)):
        group = [x for x in group if x]
        with open(f"prostate_cancer_manifest-{i}.tsv", 'w') as out:
            # out.write(head)
            for x in group:
                out.write(x)

