from Bio import SeqIO
import time
import sys
import argparse
from Bio.Align.Applications import ClustalwCommandline
import os
import re
from Bio import AlignIO
from Bio.Align import AlignInfo


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        default=sys.stdout)
    return parser


if __name__ == '__main__':
    start_time = time.time()

    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])

    pattern = re.compile("cluster_[-0-9]+.fa")
    dir_name = "."
    files = os.listdir(dir_name)
    for item in files:
        if pattern.match(item):
            cline = ClustalwCommandline("clustalw2", infile=item)
            cline()
            align = AlignIO.read(item[:-3] + ".aln", "clustal")
            summary_align = AlignInfo.SummaryInfo(align)
            consensus = summary_align.gap_consensus()
            print(item[:-3], "consensus:", consensus, file=namespace.output)
    print("Time:", time.time() - start_time)
