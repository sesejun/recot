#!/usr/bin/env python
"""Select one of the gene sequence if a gene in target species associates
with multiple genes in original species.
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
import sys, os, errno, re
import shutil
from optparse import OptionParser

def removeBoolean(id,pos,startDic):
    """
    """
    flag = False
    for k, v in startDic.items():
        dic_id = getsubString(k,'@')
        if id == dic_id:
            if pos == v:
                flag = True
                startDic.pop(k)
            else:
                continue
        #end if id == dic_id:
    #end for k, v in startDic.items():
    return flag

def removerid(sam_file, outfile, ridList,startDic):
    """
    Remove gene sequence that are overlapped with other gene sequences.
    """
    f=open(outfile, "w")
    for line in open(sam_file):
        itemList = line[:-1].split('\t')
        if itemList[0][0:1] == '@':
            f.write(line)
            continue
        id = getsubString(itemList[0], '|')
        pos = int(itemList[3])
        flag = False
        if (id in ridList) == True:
            flag = removeBoolean(id,pos,startDic)
            if flag == False:
                f.write(line)
            #end if flag == False:
        else:
            f.write(line)
    f.close()
    

def get_refFeature(target_gff_file, gene_set_dic):
    log.info('Extract gene lines in gff file')
    geneList={}#{id:start end}
    f_count = 0#feature line count
    for line in open(target_gff_file):
        itemList = line[:-1].split('\t')
        if len(itemList) != 9:
            continue
        feature = itemList[2]
        ch = itemList[0]
        '''
        if check_ch != ch:
            continue
        '''
        if feature == 'gene':
            id = itemList[8][3:itemList[8].find(';')]
            start = itemList[3]
            end = itemList[4]
            strand = itemList[6]
            if (id in gene_set_dic.values()) == True:
                geneList[id]=("\t".join([ch,start,end,strand]))
                #print id + "\t".join([start,end])
                #log.info(id+'='+geneList[id])
                f_count = f_count + 1
            else:
                continue
        elif feature != 'gene':
            continue 
        
    #end(for line in open(target_gff_file):)
    #print 'length gene list:'+ str(len(geneList))
    
    
    
    return geneList
#end (def get_refgene(target_gff_file))    

def getRemoveIds_geneSet(chrList, extend_length,gene_set_dic, working_dir, target_gff_file):
    """
    Remove gene sequence that are overlapped with other gene sequences.
    """
    sameid_c = 0#same id count
    allstartDic = {}
    ridList = []#remove id list
    file = os.path.join(working_dir, 'removeid.tab')
    test = 0
    
    refinfo_dic = get_refFeature(target_gff_file,gene_set_dic)
    
    for ch in chrList:
        startDic = {}#{id,start}
        endDic = {}#{id,end}
        sortidList = []#id
        scount = 0#stop counter
        for line in open(os.path.join(working_dir, 'chrfile.'+ ch + '.tab')):
            #print line
            itemList = line[:-1].split('\t')
            id = itemList[1]
            start = int(itemList[2])
            end = int(itemList[3])
            startDic[id+'@'+str(sameid_c)] = start
            endDic[id+'@'+str(sameid_c)] = end
            sameid_c = sameid_c + 1
            
        #sort number by start positions
        count = 0
        # we might use argsort
        sortidList = [k for (k,v) in sorted(startDic.items(), key=lambda x:x[1])]
        #for (k,v) in sorted(startDic.items(), key=lambda x:x[1])
        #    #print k,v
        #    sortidList.append(k)
        #    count = count + 1
        #end(for k,v in sorted(startDic.items(), key=lambda x:x[1]):
        
        
        
        
        #compare id
        x = 0
        y = 0
        original_id = ''
        last_idx = len(sortidList)
        while x < last_idx:
            remove_flag = 0# not remove:0 , remove 1
            tid = sortidList[x]#original id
            tstart = startDic[tid]#original start
            tend = endDic[tid]#original end
            max_region = 0
            max_regionID = ''
            max_start = 0
            max_end = 0
            y = x + 1
            
            
            # if gene_set_file has original id,
            # we have to change flag(flag name: remove_flag).
            #<flag>
            #1: This gene_set_file has id in gene sets
            #2: This gene_set_file dosen't have id in gene sets
            
            original_id = tid[:tid.find('@')]
            if (tid[:tid.find('@')] in gene_set_dic) == True:
                if (gene_set_dic[original_id] in refinfo_dic) == True:
                    log.info('gene set! (' + original_id + ', ' + gene_set_dic[original_id] + ')')
                    ref_itemList = refinfo_dic[gene_set_dic[original_id]][:-1].split('\t')
                    #log.info('original:(' + ch + ',' + str(tstart) + ',' + str(tend) + ')->ref:(' + ref_itemList[0] + ',' + ref_itemList[1] + ',' + ref_itemList[2] + ')') 
                    if ch == ref_itemList[0]:
                        refstart = int(ref_itemList[1])
                        refend = int(ref_itemList[2])
                        refstrand = ref_itemList[3]
                        if tstart > refend :
                            next
                        elif tend < refstart:
                            next
                        else:
                            remove_flag = 1
                    
                    #print 'remove_flag = ' + str(remove_flag)
                    
                #elif (gene_set_dic[original_id] in refinfo_dic) == False:
                    #log.info('tid:'+ tid + 'Original ID:' + original_id + '(' + str(tstart) + ',' + str(tend) + ')')
                    #log.info('Gene set! ...But no gff line.' + gene_set_dic[original_id] )
            elif (tid[:tid.find('@')] in gene_set_dic) == False:
                #log.info('tid:'+ tid + 'Original ID:' + original_id + '(' + str(tstart) + ',' + str(tend) + ')')
                log.info('No gene set!')
                
                
                
            while y < last_idx:
                
                cid = sortidList[y]
                cstart = startDic[cid]
                cend = endDic[cid]
                overlap_len = tend - cstart
                if overlap_len < 0:
                    # no overlapped region between tid and cid
                    break
                
                
                if tstart < cstart and cend < tend:
                    # if cid is completely inside of tid,
                    # remove tid because it may have long intron.
                    # However, this procedure might cause the problem
                    # when it has very short mapped region.
                    # We have to change the algorithm to select the best one
                    
                    
                    
                    
                   
                    if remove_flag == 1:
                        if (getsubString(cid,'@') in ridList) == False:
                            allstartDic[cid] = tstart
                            ridList.append(getsubString(cid,'@'))
                    
                    elif remove_flag == 0:
                        if (getsubString(tid,'@') in ridList) == False:
                            allstartDic[tid] = tstart
                            ridList.append(getsubString(tid,'@'))
                            #log.info('remove:' + getsubString(tid,'@'))
                    
            
                elif (overlap_len > 2*extend_length + MARGIN) or (float(overlap_len)/float(tend-tstart) > 0.5):
                    # tail of "tid" is overlapped with head of cid
                    
                    if remove_flag == 1:
                        if (getsubString(cid,'@') in ridList) == False:
                            allstartDic[cid] = cstart
                            ridList.append(getsubString(cid,'@'))
                        
                    elif remove_flag == 0:
                        if (getsubString(cid,'@') in ridList) == False:
                            allstartDic[cid] = cstart
                            ridList.append(getsubString(cid,'@'))
                            #log.info('remove:' + getsubString(cid,'@'))
                elif tend < cstart:
                    break
                
                y += 1
            
            x += 1
        ridList = list(set(ridList))
        #break
    #end(for ch in chrList:)
    
    f=open(file, "w")
    for x in ridList:
        f.write(x)
        f.write('\n')
    f.close()
    return ridList,allstartDic

def getRemoveIds_no_geneSet(chrList, extend_length, working_dir):
    """
    Remove gene sequence that are overlapped with other gene sequences.
    """
    sameid_c = 0#same id count
    allstartDic = {}
    ridList = []#remove id list
    file = os.path.join(working_dir, 'removeid.tab')
    test = 0
    for ch in chrList:
        startDic = {}#{id,start}
        endDic = {}#{id,end}
        sortidList = []#id
        scount = 0#stop counter
        for line in open(os.path.join(working_dir, 'chrfile.'+ ch + '.tab')):
            #print line
            itemList = line[:-1].split('\t')
            id = itemList[1]
            start = int(itemList[2])
            end = int(itemList[3])
            startDic[id+'@'+str(sameid_c)] = start
            endDic[id+'@'+str(sameid_c)] = end
            sameid_c = sameid_c + 1
            
        #sort number by start positions
        count = 0
        # we might use argsort
        sortidList = [k for (k,v) in sorted(startDic.items(), key=lambda x:x[1])]
        #for (k,v) in sorted(startDic.items(), key=lambda x:x[1])
        #    #print k,v
        #    sortidList.append(k)
        #    count = count + 1
        #end(for k,v in sorted(startDic.items(), key=lambda x:x[1]):
            
        #compare id
        x = 0
        y = 0
        last_idx = len(sortidList)
        while x < last_idx:
            tid = sortidList[x]#original id
            tstart = startDic[tid]#original start
            tend = endDic[tid]#original end
            y = x + 1
            while y < last_idx:
                cid = sortidList[y]
                cstart = startDic[cid]
                cend = endDic[cid]
                overlap_len = tend - cstart
                if overlap_len < 0:
                    # no overlapped region between tid and cid
                    break
                if tstart < cstart and cend < tend:
                    # if cid is completely inside of tid,
                    # remove tid because it may have long intron.
                    # However, this procedure might cause the problem
                    # when it has very short mapped region.
                    # We have to change the algorithm to select the best one
                    allstartDic[tid] = tstart
                    ridList.append(getsubString(tid,'@'))
                    
                elif (overlap_len > 2*extend_length + MARGIN) or (float(overlap_len)/float(tend-tstart) > 0.5):
                    # tail of "tid" is overlapped with head of cid
                    allstartDic[cid] = cstart
                    ridList.append(getsubString(cid,'@'))
                elif tend < cstart:
                    break
                y += 1
            x += 1
        ridList = list(set(ridList))
    #end(for ch in chrList:)
    
    f=open(file, "w")
    for x in ridList:
        f.write(x)
        f.write('\n')
    f.close()
    return ridList,allstartDic
    

def getsubString(w, c):
    """
    Extract a subequence before character 'c' appeared.
    """
    count = 0
    for x in w:
        #print x
        if x == c:
            break
        count=count+1
    return w[:count]
#def getsubString(w, c):

def makeChrFile(chrList, file, working_dir):
    """
    Read gene positions
    """
    startDic = {}
    endDic = {}
    for ch in chrList:
        chrfile = os.path.join(working_dir, 'chrfile.'+ ch + '.tab')
        f=open(chrfile, "w") 
        for line in open(file):
            itemList = line[:-1].split('\t')
            fch = itemList[0]
            if fch == ch: 
                idname = itemList[8][3:]
                start = itemList[3]
                end = itemList[4]
                f.write("\t".join([ch, idname, start, end]) + '\n')
        f.close()
    

def samToGFF(sam_file, gff_uniq_file, target_genome):
    """
    Read SAM file, extract gene names and regions, and generate their GFF file.
    """
    f=open(gff_uniq_file, "w")
    idList = []
    fileline = 0 # for debug
    startDic = {}
    endDic = {}
    chrList = []
    for line in open(sam_file):
        fileline = fileline + 1
        if line[0] == '#':
            continue
        if line[0] == '@':
            continue
        itemList = line[:-1].split('\t')
        csum = 0
        if itemList[2] == '*':
            continue
        #log.info("ID=" + itemList[0])
        ids = itemList[0].split("|")
        idname = ids[0]
        idList.append(idname)
                
        flag = itemList[1]
        rname = itemList[2]
        pos = int(itemList[3])
        mapq = itemList[4]
        cigar = itemList[5]
        seq = itemList[9]
        chrList.append(rname)
                
        a = 0
        b = 0
        seqnum = 0
        csum = pos
        idnum = 0
        count = 0
        cigpos = 0

        for x in cigar:
            op = ''
            if len(cigar) < b:
                break
            if x =='M':
                b = count
                q = int(cigar[a:b])
                idnum = idnum + 1
                seqnum = seqnum+q
                a=b+1
                b = b+1
                csum = csum+q
                #print '--------------'
            elif x == 'I':
                b = count
                q = int(cigar[a:b])
                idnum = idnum + 1
                seqnum = seqnum+q
                a=b+1
                b = b+1
                #print '--------------'
            elif x == 'D':
                b = count
                q = int(cigar[a:b])
                a=b+1
                b = b+1
                csum = csum+q
                #print '--------------'
            elif x == 'N':
                b = count
                q = int(cigar[a:b])
                a=b+1
                b = b+1
                csum = csum+q
                #print '--------------'
            elif x == 'S':
                b = count
                q = int(cigar[a:b])
                seqnum = seqnum+q
                a=b+1
                b = b+1
                #print '--------------'
            elif x == 'H':
                b = count
                q = int(cigar[a:b])
                seqnum = seqnum+q
                a=b+1
                b = b+1
                #print '--------------'
            elif x == 'P':
                b = count
                q = int(cigar[a:b])
                a=b+1
                b = b+1
                csum = csum+q
                #print '--------------'
            elif x == '=':
                b = count
                q = int(cigar[a:b])
                idnum = idnum + 1
                seqnum = seqnum+q
                a=b+1
                b = b+1
                csum = csum+q
                #print '--------------'
            elif x == 'X':
                b = count
                q = int(cigar[a:b])
                idnum = idnum + 1
                seqnum = seqnum+q
                a=b+1
                b = b+1
                csum = csum+q
                #print '--------------'
            else:
                b = b + 1 
            count = count + 1
            #break
            #print 'id=',id,  '(start, end):', pos, csum
        f.write("\t".join([rname,target_genome,'gene',str(pos), str(csum),
                           '.', '.', '.', 'ID='+idname]) + '\n')
    f.close()    
    #Compare(chrList, gff_uniq_file)
    chrList = list(set(chrList))
    chrList.sort()
    return chrList
                
def addRnameList(rnameList, itemList):
    """
    Add new chromosome to a list of target chromosomes
    """
    for x in itemList:
        if x[0:3] == 'SN:':
            # ToDo: SN is not always first place.
            # Better to use regular expressions.
            ch = x[3:]
            if (ch in rnameList) == False:
                rnameList.append(ch)
                break
    return rnameList

def getRefName(read_file,rnameList):
    """
    Extract a list of target chromosomes
    """
    lineList = []
    for line in open(read_file):
        itemList = line[:-1].split('\t')
        if itemList[0][0:1] == '@':
            if itemList[0] == '@SQ':
                rnameList = addRnameList(rnameList, itemList)
                continue
        else:
            break
    #end(for line in open(file):)
    log.info("Number of Target Chromosomes/Genes: " + str(len(rnameList)))
    return rnameList

def get_geneSetDic(gene_set_file):
    """
    Get dictionary of gene sets
    """
    dic={}#{original speciese: main speciese}
    for line in open(gene_set_file):
        #print line
        itemList = line[:-1].split('\t')
        original_id = itemList[0] #original species ID
        main_id = itemList[1] #main species ID
        if main_id != '': 
            dic[original_id] = main_id
    #end(for line in open(gene_set_file):)
    
    return dic
#end (get_geneSetDic(gene_set))


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

    global MARGIN
    MARGIN = 100

    config_file = opt.config
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    try:
        log.info("Importing settings from %s" % config_file )
        extend_length = int(config.get("extract", "extend_length"))
        log.info("extend_length=" + str(extend_length))
        target_genome = config.get("global","target_genome")
        log.info("target_genome=" + target_genome) 
        target_gff_file = config.get("global","target_gff_file")
        log.info("target_genome=" + target_gff_file)
        gene_seq_sam = config.get("genemap", "gene_seq_sam")
        log.info("gene_seq_sam=" + gene_seq_sam)
        gene_seq_uniq = config.get("removemultiple","gene_seq_uniq_sam")
        log.info("gene_seq_new" + gene_seq_uniq)
        gff_uniq_file = config.get("removemultiple", "gff_uniq_file")
        log.info("gff_uniq_file=" + gff_uniq_file)
        gene_set = config.get("removemultiple","gene_set")
        log.info("gene_set=" + gene_set)
        working_dir = config.get("removemultiple", "working_dir")
        # To show next command
        original_fasta_file = config.get("global", "original_fasta_file")
        target_fasta_file = config.get("global", "target_fasta_file")
        read_sam_file=config.get("combine_reads", "read_sam_file")
    except ConfigParser.NoOptionError:
        print "Option name missing. Check your setting.ini file"
        raise
    # make working directory if it doesn't exist.
    try:
        os.makedirs(working_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

    chrList = []
    geneset_dic = {}
    log.info("Convert SAM file to GFF")
    chrList = samToGFF(gene_seq_sam, gff_uniq_file, target_genome)
    log.info("Extract Chromosome Names")
    getRefName(gene_seq_sam, chrList)
    log.info("Extract Chromosome Names")
    makeChrFile(chrList, gff_uniq_file, working_dir)
    log.info("Find Overlapped Regions")
    if gene_set == 'no':
        log.info("Not use gene sets")
        ridList,startDic = getRemoveIds_no_geneSet(chrList, extend_length, working_dir)
    else:
        log.info("Use gene sets")
        geneset_dic = get_geneSetDic(gene_set)
        geneset_size = len(geneset_dic)
        log.info("Number of gene sets:"+str(geneset_size))
        ridList,startDic = getRemoveIds_geneSet(chrList, extend_length,geneset_dic, working_dir, target_gff_file)
        
        
    log.info("Remove Overlaps")
    removerid(gene_seq_sam, gene_seq_uniq,ridList,startDic)
    log.info("Result File is " + gene_seq_uniq)
    
    #shutil.rmtree(working_dir)

    print("\n=====")
    print("Generated new SAM format file: %s" % gene_seq_uniq)
    print("Then, you prepare files to change coordinages by running:")
    print("$ python combine_reads_genome.py")
    print("")
    print("Before the step, if you have not map your reads agains original genes/genomes,")
    print("you have to do it now.")
    print("\tEg.")
    print("$ bwa aln %s your_read.fastq > your_read.sai" % (target_fasta_file))
    print("$ bwa samse %s your_read.sai your_read.fastq > %s" % (original_fasta_file, read_sam_file))


if __name__ == "__main__": 
    #param = sys.argv
    main()
