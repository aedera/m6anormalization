import argparse
import multiprocessing
import numpy as np
import time

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

class NormalizationCalculator(multiprocessing.Process):
    def __init__(self, outname, queue):
        super().__init__()
        self.outname = outname
        self.q = queue

    def run(self):
        msums, mcounts = {}, {}
        while True:
            e = self.q.get()
            if e is None:
                self.q.put(None)
                break

            sums, counts = e
            assert len(sums) == len(counts)

            for kmer in sums:
                s, c = sums[kmer], counts[kmer]
                if kmer not in msums:
                    msums[kmer]   = 0
                    mcounts[kmer] = 0
                msums[kmer]   += s
                mcounts[kmer] += c

        fout = open(self.outname, 'w')
        for kmer in msums:
            norm = msums[kmer]/mcounts[kmer]
            fout.write(f'{kmer}\t{norm:.4f}\n')
        fout.close()

class Counter(multiprocessing.Process):
    def __init__(self, idx, queue, oqueue):
        super().__init__()
        self.idx = idx
        self.q = queue
        self.o = oqueue

    def run(self):
        print(f'Starting Counter {self.idx}')
        sums, counts   = {}, {}
        while True:
            e = self.q.get()
            if e is None:
                self.q.put(None)
                break

            for kmer, mlevel in e:
                if kmer not in sums:
                    sums[kmer]   = 0.
                    counts[kmer] = 0.
                sums[kmer]   += mlevel
                counts[kmer] += 1

            if len(sums) > 1_000:
                self.o.put([sums, counts])
                sums, counts = {}, {}

            if self.o.qsize() > 1_000:
                # print(f'Counter {self.idx} is sleeping')
                time.sleep(1)
        if len(sums) > 0:
            self.o.put([sums, counts])
        print(f'Stop Counter {self.idx}')

class BEDReader(multiprocessing.Process):
    def __init__(self, chrn, fas_file, bed_file, queue):
        super().__init__()
        self.chrn = chrn
        self.fas  = fas_file
        self.bed  = bed_file
        self.q    = queue

    # def get_kmer(self, p, wsize=6):
    #     kmer = self.seq[p-wsize:p+wsize+1]
    #     if len(kmer) < (wsize*2) + 1:
    #         str_kmer = None
    #     else:
    #         str_kmer = ''.join(kmer)
    #         focal_pos = kmer[wsize]

    #         if focal_pos != 'A' and focal_pos != 'T':
    #             str_kmer = None

    #         if focal_pos == 'T': # minus strand, calculate reverse complement
    #             str_kmer = str_kmer.translate(self.tr)[::-1]    

    #     return str_kmer

    def run(self):
        print(f'Starting Reader for {self.chrn}')

        # ForkedPdb().set_trace()

        # read target chromosome from fasta file
        self.seq = read_fas(self.fas, self.chrn)

        buff = []
        with open(self.bed) as f:
            for a in f:
                _chrn, start, _, mlevel, _, cov = a.strip().split('\t')
                # _chrn = f'Chr{_chrn}'
                if _chrn != self.chrn:
                    continue

                start  = int(start)
                cov    = int(cov)
                mlevel = float(mlevel)

                if cov <= 0: # exclude positions without coverage
                    continue

                kmer = get_kmer(self.seq, start)

                buff.append([kmer, mlevel])
                if len(buff) > 1_000:
                    # print(f'{self.chrn} adding buf')
                    self.q.put(buff)
                    buff = []

                while self.q.qsize() > 1_000:
                    # print(f'sleeping {self.chrn}')
                    time.sleep(1)
        if len(buff) > 0:
            self.q.put(buff)
            buff = []
        print(f'Stop Reader for {self.chrn}')

def main(args):
    bed_fname  = args.bed
    fas_fn     = args.fas
    chrs       = args.chrs
    outname    = args.out

    chrs = [chrn.strip() for chrn in chrs.split(',')]

    iqueue = multiprocessing.Queue() # input queue
    oqueue = multiprocessing.Queue() # output queue

    # process responsable for writing the kmer and their mean
    # methylation value
    nc = NormalizationCalculator(outname, oqueue)
    nc.start()

    num_counters = len(chrs) * args.num_counters
    counters = []
    for i in range(NUM_COUNTER):
        counter = Counter(i, iqueue, oqueue)
        counter.start()
        counters.append(counter)



    processes = []
    for chrn in chrs: # launch a process for each chromosome
        p = BEDReader(chrn, fas_fn, bed_fname, iqueue)
        p.start()
        processes.append(p)

    [p.join() for p in processes]
    iqueue.put(None)
    [p.join() for p in counters]
    oqueue.put(None)
    nc.join()

def argparser():
    parser = argparse.ArgumentParser(
        'Genereate normalization constants for kmers',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument(
        "--bed",  
        required=True, 
        type=str, 
        help='bed file containing m6a levels called from Nanopore reads'
    )
    parser.add_argument(
        "--fas",  
        required=True, 
        type=str, 
        help='fasta file containing the reference genome used for m6a calling'
    )
    parser.add_argument(
        "--chrs", 
        required=True, 
        type=str, 
        help='comma-separated list of the chromosomes used for normalization'
    )
    parser.add_argument(
        "--out",  
        required=True, 
        type=str, 
        help='output file containing kmer normalization constants'
    )
    parser.add_argument(
        "--num_counters",  
        required=False, 
        type=int, 
        default=10, 
        help='number of processes used for processing each chromosome'
    )
    return parser
