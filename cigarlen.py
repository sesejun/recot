#!/usr/bin/env python
import re, sys

def cigarLength(cigar):
    p = re.compile("([0-9]+)([MIDNSHP=X])")
    len = 0
    n_of_nuc = 0
    next_cigar = cigar
    m = p.match(next_cigar)
    while m:
        n_of_nuc = int(m.group(1))
        operation = m.group(2)
        #print str(n_of_nuc) + ":" + operation
        if operation in ["M","I","S", "=", "X"]:
            len += n_of_nuc
        next_cigar = next_cigar[m.end():]
        m = p.match(next_cigar)
    return len

for line in sys.stdin:
    cols = line[:-1].split("\t")
    if len(cols) < 9:
        continue
    print str(cigarLength(cols[5])) + "\t" + cols[5] + "\t" + cols[0]
