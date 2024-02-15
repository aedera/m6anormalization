import argparse
import multiprocessing

from m6anormalization.utils.io import read_fas
from m6anormalization.utils.utils import get_kmer

# import pdb

# class ForkedPdb(pdb.Pdb):
#     """A Pdb subclass that may be used
#     from a forked multiprocessing child

#     """
#     def interaction(self, *args, **kwargs):
#         _stdin = sys.stdin
#         try:
#             sys.stdin = open('/dev/stdin')
#             pdb.Pdb.interaction(self, *args, **kwargs)
#         finally:
#             sys.stdin = _stdin

class Writer(multiprocessing.Process):
    def __init__(self, knorm_fn, fas_fn, bed_fn, chrn):
        super().__init__()
        self.knorm_fn = knorm_fn
        self.fas_fn   = fas_fn
        self.bed_fn   = bed_fn
        self.tchrn    = chrn

    def run(self):
        # read kmer normalization
        k2v = {}
        with open(self.knorm_fn) as f:
            for line in f:
                kmer, val = line.strip().split('\t')
                k2v[kmer] = float(val)

        # read fasta file to retrieve k-mers associated to positions
        seq = read_fas(self.fas_fn, self.tchrn)

        # process input bed file
        out = open(f'{self.tchrn}.bed', 'w')
        # process input bed file
        with open(self.bed_fn) as f:
            for line in f:
                chrn, start, stop, mlevel, strand, cov = line.strip().split('\t')
                if chrn != self.tchrn:
                    continue
                start  = int(start)
                mlevel = float(mlevel)

                kmer = get_kmer(seq, start)
                norm = k2v[kmer] if kmer in k2v else None

                if not norm is None:
                    if norm != 0:
                        new_level = mlevel/norm
                    else:
                        new_level = 0.

                    out.write(f'{self.tchrn}\t{start}\t{start+1}\t{mlevel}\t{strand}\t{cov}\t{kmer}\t{new_level:.4f}\n')
        out.close()

def main(args):
    knorm_fn = args.norm
    fas_fn   = args.fas
    bed_fn   = args.bed
    chrs     = args.chrs
    chrs     = [chrn.strip() for chrn in chrs.split(',')]

    processes = []
    for chrn in chrs: # launch a process for each chromosome
        p = Writer(knorm_fn, fas_fn, bed_fn, chrn)
        p.start()
        processes.append(p)
    [p.join() for p in processes]

def argparser():
    parser = argparse.ArgumentParser(
        'Apply normalization constants to genome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("--norm", required=True, type=str)
    parser.add_argument("--fas",  required=True, type=str)
    parser.add_argument("--bed",  required=True, type=str)
    parser.add_argument("--chrs", required=True, type=str)
    return parser
