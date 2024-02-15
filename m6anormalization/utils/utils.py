TR = str.maketrans('ACGTN', 'TGCAN')

def get_kmer(seq, p, wsize=6):
    kmer = seq[p-wsize:p+wsize+1]
    if len(kmer) < (wsize*2) + 1:
        str_kmer = None
    else:
        str_kmer = ''.join(kmer)
        focal_pos = kmer[wsize]

        if focal_pos != 'A' and focal_pos != 'T':
            str_kmer = None

        if focal_pos == 'T': # minus strand, calculate reverse complement
            str_kmer = str_kmer.translate(TR)[::-1]    

    return str_kmer
