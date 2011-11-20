#!/usr/bin/env python
"""Split reads and annotaions into each chromosome/gene/contig file.
This is a preparation to change coordinate of reads.
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

import os, re, sys, logging
import ConfigParser
import shutil
import errno
from optparse import OptionParser


def getsubString(w, c):
    """
    Extract string before character c from w
    """
    count = 0
    for x in w:
        if x == c:
            break
        count=count+1
    return w[:count]
#def getsubString(w, c):

def makeGenefiles(gff_file,geneList,rnameList):
    '''
    Split gene annotations in GFF file into each chromome file
    '''
    size = len(rnameList)
    prev = 0
    ends = range(0, size, 20)
    ends += [size]
    ends.pop(0)
    p = re.compile("ID=([^;]+);")
    
    f_gff = open(gff_file)
    data1 = f_gff.read()
    f_gff.close()
    
    lineLists = []
    for line in data1.split('\n'):
        itemList = line[:-1].split('\t')
        if len(itemList) < 9:
            continue
        if itemList[2] != 'gene':
            continue
        lineLists.append(itemList)
    
    for i in ends:
        chrs = rnameList[prev:i]
        f = []
        for j in range(0,i-prev):
            genefile = os.path.join(working_dir, 'removeOverlap.'+chrs[j]+'.gff')
            log.info('Generating ' + genefile)
            f.append(open(genefile, "w"))
        for itemList in lineLists:
            for j in range(0,i-prev):
                if chrs[j] == itemList[0]:
                    m = p.match(itemList[8])
                    id = ""
                    if m:
                        id = m.group(1)
                    else:
                        id = itemList[8]
                    if (id in geneList) == True:
                        f[j].write("\t".join(itemList) + "\n")
        for fp in f:
            fp.close()
        prev = i


def getGeneId(gene_sam_file):
    '''
    Extract gene name list from sam file
    '''
    rnameList = []
    geneList = []
    for line in open(gene_sam_file):
        itemList = line[:-1].split('\t')
        if itemList[0] == '@SQ':
            continue
        id = getsubString(itemList[0],'|')
        #flag = itemList[1]
        #rname = itemList[2]
        #pos = int(itemList[3])
        #mapq = itemList[4]
        #cigar = itemList[5]
        #seq = itemList[9]
        if (id in geneList) == False:
            geneList.append(id)
    return geneList

def getrnaList(rnameList, itemList):
    '''
    Parse chromosome name
    '''
    for x in itemList:
        if x[0:3] == 'SN:':
            ch = x[3:]
            if (ch in rnameList) == False:
                rnameList.append(ch)
                break
    return rnameList


def getRefName(read_file):
    '''
    Extract chromosome names
    '''
    rnameList = []
    log.info("Extract Target Chromosome Name")
    for line in open(read_file):
        itemList = line[:-1].split('\t')
        if itemList[0][0:1] == '@':
            if itemList[0] == '@SQ':
                rnameList = getrnaList(rnameList, itemList)
                continue
        else:
            break
    log.info("# of target seqs: " + str(len(rnameList)))
    return rnameList

def getReadSamFile(read_file,rnameList):
    """
    Separate sam file into each chromosome file
    """
    size = len(rnameList)
    prev = 0
    ends = range(0, size, 20)
    ends += [size]
    ends.pop(0)
    
    
    
    for i in ends:
        chrs = rnameList[prev:i]
        f = []
        ch_p = ''
        jj = 0
        for j in range(0,i-prev):
            samfile = os.path.join(working_dir, 'MappedRead.'+chrs[j]+'.sam')
            log.info('Generating ' + samfile)
            f.append(open(samfile, "w"))
        for line in open(read_file, "r"):
            
            itemList = line[:-1].split('\t')
            
            if len(itemList) < 11:
                continue
            #print itemList
            if itemList[0][0:1] == '@':
                continue
            line_ch = itemList[2]
            if line_ch == '*':
                continue
            if int(itemList[1]) & 0b100 != 0:
                continue
            
            if ch_p != line_ch:
                for j in range(0,i-prev):
                    if chrs[j] == line_ch:
                        f[j].write(line)
                        jj = j
                        ch_p = line_ch
                        continue
                #end for j in range(0,i-prev):
            elif ch_p == line_ch:
                f[jj].write(line)
            '''
            for j in range(0,i-prev):
                if chrs[j] == line_ch:
                    f[j].write(line)
                    continue
            '''
        for fp in f:
            fp.close()
        prev = i


def getGFFStartEnd(file, len_param):
    """
    Extract start and end positions of genes
    """
    dicS = {}
    dicE = {}
    direct = {}
    for line in open(file):
        itemList = line[:-1].split('\t')
        start = int(itemList[3])-len_param
        if start <0:
            start = 0
        end = int(itemList[4])+len_param
        #id = getsubString(itemList[8][4:],';') # ToDo: need to check for other species
        id = itemList[8][itemList[8].find('=')+1:itemList[8].find(';')]
        dicS[id]= start
        dicE[id]= end
        direct[id] = itemList[6]
    return dicS,dicE,direct

def sortId(dicS):
    """
    Sorting id according to positions
    """
    idList = {}
    count = 0
    for k,v in sorted(dicS.items(), key=lambda x:x[1]):
        idList[count]=k
        count = count + 1
    return idList


def getSeqLength(cigar):
    '''
    Calculate sequence length from cigar sequence
    '''
    length = 0
    
    '''
    CIGAR====================================================
    M 0 alignment match (can be a sequence match or mismatch) 
    I 1 insertion to the target 
    D 2 deletion from the target 
    N 3 skipped region from the target 
    S 4 soft clipping (clipped sequences present in SEQ) 
    H 5 hard clipping (clipped sequences NOT present in SEQ) 
    P 6 padding (silent deletion from padded target) 
    = 7 sequence match 
    X 8 sequence mismatch 
    ========================================================

    '''
    a = 0
    b = 0
    seqnum = 0
    idnum = 0
    count = 0
    cigpos = 0
    cigar_len = len(cigar)
    
    for x in cigar:
        if cigar_len < b:
            break
        if x >= '0' and x <= '9':
            b = b + 1
        else:
            b = count
            
            q = int(cigar[a:b])
            if x =='M':
                idnum = idnum + 1
                seqnum = seqnum + q
                length = length + q
            elif x == 'D' or x == 'N' or x == 'S' or x == 'H':
                length = length+q
            elif x == 'I' or x == '=' or x == 'X':
                idnum = idnum + 1
                seqnum = seqnum+q
            #elif x == 'P':
                # Nothing to do
            a = b+1
            b = b+1
        count = count + 1
    return length

def getSAMStartEnd(file):
    """
    
    """
    readS = {}
    readE = {}
    readDic = {}
    for line in open(file):
        #print line
        itemList = line[:-1].split('\t')
        id = itemList[0]
        flag = int(itemList[1])
        start = int(itemList[3])
        cigar = itemList[5]
        #seq = itemList[9]
        if flag & 0b100 != 0:
            continue
        seq_len = getSeqLength(cigar)
        readS[id] = start
        readE[id] = start + seq_len -1
        readDic[id] = line
    return readS,readE,readDic

def getsameIDList(id, file):
    """
    Read Gene ID List
    """
    glineList = []
    newread = []
    
    for line in open(file):
        itemList = line[:-1].split('\t')
        line_id = getsubString(itemList[0],'|')
        
        if id == line_id:
            glineList.append(line)
        else:
            newread.append(line)
        return glineList

def getReadOnGeneFile(rnameList, len_param):
    """
    Select reads that are on gene
    """
    log.info("Select reads that are on genes")
    for ch in rnameList:
        tcount = 0
        
        geneS = {}#gene start
        geneE = {}#gene end
        g_direct = {}#gene direction
        readS = {}#read start
        readE = {}#read End
        readDic = {}#readDic[id] = read
        sortGeneId = {}
        sortReadId = {}
        genefile = os.path.join(working_dir, 'removeOverlap.'+ch+'.gff')
        readfile = os.path.join(working_dir, 'MappedRead.'+ch+'.sam')
        rgfile = os.path.join(working_dir, 'ReadOnGeneList.'+ch+'.tab')
        log.info("Generate " + rgfile)
        f=open(rgfile, "w") 
        
        geneS, geneE, g_direct = getGFFStartEnd(genefile, len_param)
        sortGeneId = sortId(geneS)
        
        readS, readE,readDic = getSAMStartEnd(readfile)
        sortReadId = sortId(readS)
        ys = 0
        
        for x in range(len(sortGeneId)):
            
            gID = sortGeneId[x]#gene id
            gs = geneS.get(gID)#gene start
            ge = geneE.get(gID)#gene end
            gd = g_direct.get(gID)
            glineList = []
            sameG = False
            
            for y in range(ys,len(sortReadId)):
                rID = sortReadId[y]
                rs = readS.get(rID)
                re = readE.get(rID)
                if rs >= gs:
                    if re <= ge:
                        f.write(gID)
                        f.write('\t')
                        f.write(str(gs))
                        f.write('\t')
                        f.write(str(ge))
                        f.write('\t')
                        f.write(gd)
                        f.write('\t')
                        f.write(rID)
                        f.write('\t')
                        f.write(str(rs))
                        f.write('\t')
                        f.write(str(re))
                        f.write('\t')
                        f.write(readDic.get(rID))
                    elif re > ge:
                        ys = y
                        break
                elif rs > ge:
                    ys = y
                    break
        f.close()
    
def main():
    """
    Main Function
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
    logging.basicConfig(level=logging.DEBUG,
                        format="%(name)s: %(levelname)s @ %(asctime)s %(message)s")

    CONFIG_FILE = 'settings.ini'

    config_file = opt.config
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    global working_dir
    try:
        log.info("Importing settings from %s" % config_file )
        gene_sam_file = config.get("removemultiple", "gene_seq_uniq_sam")
        log.info("gene_sam_file=" + gene_sam_file)
        gff_file = config.get("extract", "gff_file")
        log.info("gff_file=" + gff_file)
        read_sam_file = config.get("combine_reads", "read_sam_file")
        log.info("read_sam_file=" + read_sam_file)
        #final_sam_file = config.get("combine_reads", "final_sam_file")
        #log.info("final_sam_file=" + final_sam_file)
        working_dir = config.get("combine_reads", "working_dir")
        log.info("working_dir=" + working_dir)
        len_param = int(config.get("extract", "extend_length"))
        log.info("len_param=" + str(len_param))
    except ConfigParser.NoOptionError:
        print "Option name missing. Check your setting.ini file"
        raise

    try:
        os.makedirs(working_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

    geneList = []
    rnameList = [] # Target Chromosome Names
    
    rnameList = getRefName(read_sam_file)
    getReadSamFile(read_sam_file,rnameList)
    
    geneList = getGeneId(gene_sam_file)

    makeGenefiles(gff_file,geneList,rnameList)
    getReadOnGeneFile(rnameList,len_param)
    

    print("=====")
    print("Successfully finished.")
    print("You next convert the read coordinate by running:")
    print("$ python convert_coordinate.py")

if __name__ == "__main__": 
    main()
