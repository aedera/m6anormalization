import numpy as np

def read_fas(fname, target_chrn):
    seqs = []
    with open(fname) as f:
        for line in f:
            line = line.rstrip()
            if line[0] == '>':
                chrn = line[1:]
            else:
                if chrn == target_chrn:
                    seqs.append(line.upper())

    seqs = np.array(list(''.join(seqs)))

    return seqs
