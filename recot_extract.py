#!/usr/bin/env python
"""Extract genes and their reguratory region sequences.
Configurations have to be described in ``settings.ini''

"""
"""
Copyright (c) 2011, RECOT Development Team
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of the Database Center for Life Science, the Universiry of Tokyo, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

__author__ =  'Akiko Izawa and Jun Sese'
__version__=  '0.5'

import ConfigParser
import logging
import sys, os
import re
from optparse import OptionParser

def getChrSeq(fasta_file):
    """
    get entire DNA sequences from 'fasta_file'
    """
    seq = ""
    chseq = {}
    c = 0
    f = open(fasta_file, 'r')
    lines = f.readlines()
    p = re.compile('>(.+)')
    for line in lines:
        line = line.rstrip()
        line = line.lstrip()
        if line[0] == "#":
            next
        m = p.match(line)
        if m:
            if not seq == "":
                chseq[ch] = seq
                log.info("Chromosome length of " + ch + ": " + str(len(seq)))
                seq = ""
            ch = m.group(1)
            log.info("Reading chromosome name: " + ch)
        else:
            seq = seq + line
    if not seq == "":
        chseq[ch.split()[0]] = seq
        log.info("Chromosome length of " + ch + ": " + str(len(seq)))
    return chseq

def getID(info):
    """
    extract id from GFF information
    """
    idnum = ""
    p = re.compile('ID=([^;]+);')
    m = p.match(info)
    if m:
        idnum = m.group(1)
    else:
        log.error("No ID found in " + info)
    return idnum

def changeStrand(seq):
    """
    generate reverse strand sequence
    """
    revdir = seq[::-1]
    relation = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G',
        'a':'t',
        't':'a',
        'g':'c',
        'c':'g'
    }
    revstr = ""
    for c in revdir:
        if relation.has_key(c):
            revstr += relation[c]
        else:
            revstr += c
    return revstr


def makeGeneGFF(extlen, fasta_file,gff_file, outfile):
    """
    extract gene region sequences. The genes are identified from GFF file.
    """
    # initialize variables
    ch = ""
    id = ""
    seq = ""
    # read all chromosome
    chrseq = getChrSeq(fasta_file)
    # parse files
    f = open(outfile, "w+")
    for line in open(gff_file): 
        itemList = line[:-1].split('\t')
        # each line contains 9 items
        if len(itemList) != 9:
            continue
        else:
            # select gene region
            if itemList[2]=="gene":
                # get a gene name and its position
                ch = itemList[0]
                pre = int(itemList[3])
                pos = int(itemList[4])
                strand = itemList[6]
                info = itemList[8]
                idnum = getID(info)
                log.debug("Processing gene " + idnum + " at " + str(pre) + "," + str(pos) + " on " + ch)
                # extract gene sequence
                if (ch in chrseq) == False:
                    continue
                seq = chrseq[ch]
                startpos = max(min(pre,pos)-extlen, 0)
                endpos = min(max(pre,pos)+extlen, len(seq)-1)
                subseq = seq[startpos:endpos]
                # change strand if it's in reverse strand
                if strand == '-':
                    subseq = changeStrand(subseq)
                log.debug("Length: " + str(len(subseq)))
                f.write(">"+idnum+"|"+ch)
                f.write('\n')
                f.write(subseq)
                f.write('\n')
    f.close()

def main():
    """
    Main function
    """

    parser = OptionParser()
    parser.add_option("-d", "--debug", action="store_true", dest="debug",
                      default=False,
                      help="print debug messages")
    parser.add_option("-q", "--quiet", action="store_true", dest="quiet",
                      default=False,
                      help="quiet mode. only report errors and important messages.")
    parser.add_option("-c", "--config", action="store", dest="config",
                      default="settings.ini",
                      help="change config file [default: %default]")
    (opt, args) = parser.parse_args()


    # set up logging info.
    global log
    debug_level = logging.INFO
    if opt.debug and opt.quiet:
        print "--debug and --quiet are mutual option. Select one of them."
        return
    if opt.debug:
        debug_level = logging.DEBUG
    if opt.quiet:
        debug_level = logging.ERROR
    
    script_name = sys.argv[0].split('/')[-1]
    log = logging.getLogger(script_name)
    logging.basicConfig(level=debug_level,
                        format="%(name)s: %(levelname)s @ %(asctime)s %(message)s")
    config_file = opt.config
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    try:
        target_genome = config.get("global", "target_genome")
        extlen = int(config.get("extract", "extend_length"))
        log.info("extend_length=" + str(extlen))
        fasta_file = config.get("global", "original_fasta_file")
        log.info("fasta_file=" + fasta_file)
        gff_file = config.get("extract", "gff_file")
        log.info("gff_file=" + gff_file)
        gene_seq_file = config.get("extract", "gene_seq_file")
        log.info("gene_seq_file=" + gene_seq_file)
        gene_sam_file = config.get("genemap", "gene_seq_sam")
    except ConfigParser.NoOptionError:
        print "Option name missing. Check your setting.ini file"
        raise
    makeGeneGFF(extlen, fasta_file, gff_file, gene_seq_file)

    print("\n=====")
    print("Generated fasta format gene sequence file in %s" % gene_seq_file)
    print("with %s-bp upstream and downstream regions" % extlen)
    print("Next step is to map the sequences onto target genome.")
    print("Eg. gmap_build -d %s sample/dmel-2R.fasta" % target_genome)
    print("    gmap -f samse -d %s %s > %s" % (target_genome, gene_seq_file, gene_sam_file))
    print("Then, you remove redundant associations of genes by running:")
    if config_file == "settings.ini":
        print("$ python recot_rm_overlap.py")
    else:
        print("$ python recot_rm_overlap.py -c " + config_file)
    print("")
#end(def main():)


if __name__ == "__main__": 
    main()
