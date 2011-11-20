#!/usr/bin/env python

"""Convert coordinate of genes to target genome.
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

import shutil
import re, sys, errno, os
import logging
import ConfigParser
from optparse import OptionParser


def getChrSAMfile(gene_file, rnameList):
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
            sam_file = os.path.join(working_dir, 'GeneOnMainRef.'+chrs[j]+'.sam')
            log.info('Generating ' + sam_file)
            f.append(open(sam_file, "w"))
        for line in open(gene_file):
            itemList = line[:-1].split('\t')
            # Check: the below line is id = itemList[0].split("|")[0] and [1]?
            # or p = re.compile("([^|]+)|(.+)"); m.match(itemList[0]);
            # m.group(1) and (2) ?
            #id = getsubString_a(itemList[0],'|')
            
            if len(itemList) < 11:
                continue
            
            line_ch = getsubString_b('|',itemList[0])
                    
            
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
                
                
        for fp in f:
            fp.close()
        prev = i

def getCoorFlag(c):
    flag =  0
    if c == 'I':
        flag = 1
    elif c == 'P':
        flag = 1
    return flag
#end def getCoorFlag(c):

def getCoorNum(num,flag,coor_num,coor,count):
    if flag == 0:
        for i in range(count,count+num):
            coor[coor_num] = str(i)
            coor_num = coor_num + 1
        count = count+num
    elif flag == 1:
        for i in range(0,num):
            coor[coor_num] = '*'
            coor_num = coor_num + 1
    return count,coor_num
#end getCoor(coorflag,coorS,coorE,coor)

def getCigarDictionary(cigar,pos,num,dic):
    for i in range(0,num):
        dic[pos] = cigar
        pos=pos+1
    return pos,dic
#end(def getCigarDictionary(cigar,pos,num,dic):)

def output(seqdic,seq_num,num,c,s,e,seq,flag):#addCigarSeq
    #print num,c
    #print '(s,e)->',s,e
    if flag == 0:
        for pos in range(s,e):
            seqdic[seq_num] = seq[pos:pos+1]
            seq_num = seq_num + 1
        #end for pos in range(s,e):
    elif flag == 1:
        s = ''
        for i in range(0,int(num)):
            s = s + '*'
            seqdic[seq_num] = '*'
            seq_num = seq_num + 1
        #print s
    elif flag == 2:
        s = ''
        for i in range(0,int(num)):
            s = s + '.'
            seqdic[seq_num] = '.'
            seq_num = seq_num + 1
       # print s
    return seqdic,seq_num
#end def output(num,c,s,e,seq):

# ToDo: Hope to change function name
def getsubString_a(w, c):
    #print '* getsubString!'
    count = 0
    for x in w:
        if x == c:
            break
        count=count+1
    #end(for x in w:)
    return w[:count]
#def getsubString(w, c):

# ToDo: Hope to change function name
def getsubString_b(c,w):
    count = 0
    for x in w:
        #print x
        if x == c:
            break
        count=count+1
    #end(for x in w:)
    return w[count+1:]
#end(def getsubString('|',itemList[0]))

def getSamIDList(ch,rg_gID, gene_file):
    geneSamList = []
    #newread = []
    for line in open(gene_file):
        itemList = line[:-1].split('\t')
        line_id = getsubString_a(itemList[0],'|')
        line_ch = getsubString_b('|',itemList[2])
        #if (line_ch == '*') or (line_ch == ch):
        if line_ch == '*':
            continue
        if rg_gID == line_id:
            geneSamList.append(line)
        #else:
        #    newread.append(line)
    #end(for line in open(gene_file):)
    return geneSamList
#end(getSamIDList(id, gene_file))

def getrnaList(rnameList, itemList):
    for x in itemList:
        if x[0:3] == 'SN:':
            ch = x[3:]
            if (ch in rnameList) == False:
                rnameList.append(ch)
                break
            #end(if (ch in rnameList) == False:)
        #end(if x[0:3] == 'SN:':)
    #end(for x in itemList:)
    return rnameList
#end(def getrnaList(rnameList, itemList):)

def getRefName(read_file,rnameList):
    lineList = []
    log.info('get Target Chromosome Names')
    for line in open(read_file):
        itemList = line[:-1].split('\t')
        if itemList[0][0:1] == '@':
            if itemList[0] == '@SQ':
                rnameList = getrnaList(rnameList, itemList)
                continue
        else:
            break
    return rnameList

def getSAMheader(gene_file):
    list = []
    for line in open(gene_file):
        itemList = line[:-1].split('\t')
        if itemList[0][0:1]=='@':
            list.append(line)
    return list

def reverseBoolean(flag):
    int_flag = int(flag)
    boolean = False
    if int_flag & 0b10000 != 0:
        boolean = True
    return boolean
#end def reverseBoolean():

def addCigarSeq(cigar,seq):
    num = ''#cigar nomber(String)
    int_num = 0 #cigar number(int)
    s = 0# seq start position
    e = 0# seq end position
    
    coordic = {}#(String)
    subcoor = ''
    coorS = 0#coor start(int)
    coorE = 0#coor end(int)
    coor_c = 0 #coor count(int)
    coorflag = 0#0->num,1->'*' 
    coordic_num = 0
    
    cigardic = {}
    cigar_num = 0
    
    seqdic = {}
    seq_num = 0
    
    for c in cigar:
        if c >= '0' and c <= '9':
            num = num + c
            continue
        int_num = int(num)
        coorflag = getCoorFlag(c)
        if c == 'M':
            e = int_num + s
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,0)
        elif c == 'I':
            e = int_num + s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,0)
        elif c == 'D':
            e = s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,1)
        elif c == 'N':
            e = s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,2)
        elif c == 'S':
            e = int_num + s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,0)
        elif c == 'H':
            e = int_num + s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,1)
        elif c == 'P':
            e = s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,1)
        elif c == '=':
            e = int_num + s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,0)
        elif c == 'X':
            e = int_num + s
            coorS = s
            coorE = e
            coor_c,coordic_num = getCoorNum(int_num,coorflag,coordic_num,coordic,coor_c)
            cigar_num,cigardic = getCigarDictionary(c,cigar_num,int_num,cigardic)
            seqdic,seq_num = output(seqdic,seq_num,num,c,s,e,seq,0)
        s = e
        num = ''
             
        
    #end for c in cigar:
    
    
        
    return coordic,cigardic,seqdic
#end def addCigarSeq(cigar,seq):

def Get_gene_spos(g_pos,gs,rs,gene_coorDic,gene_cigarDic,gene_seqDic):
    diff = rs - gs
    basecigar = ''
    newcigar = ''
    newpos = 0
    gene_spos = 0
    seq_spos = 0
    
    coorflag = 0#0:number,1:'*'
    
    seq_num= 0
    for gc in range(0,len(gene_coorDic)):
        
        
            
        
        if (gene_seqDic[gc] == '.') or (gene_seqDic[gc] == '*'):
            #print 'seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]',seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc],'----skip--------'
            continue
        else:
            #print 'seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]',seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc],'gene_spos:',gene_spos
            if seq_spos == diff:
                break
            if gc == 0:
                gene_spos = gc
            else:
                gene_spos = gc
            
            
                
                
                
            
                    
            #end if seq_spos == diff:
            seq_num = gc
            seq_spos = seq_spos + 1
            #end if gene_spos == diff:
        # end else
        
    #end for gc in range(diff,len(gene_coorDic)):
    start_base = 0
    for gc in range(0,len(gene_coorDic)):
        if (gene_coorDic[gc]=='*') or (gene_cigarDic[gc]=='S'):
            
            continue
        else:
            newpos = g_pos+start_base
            start_base = start_base + 1
            if gene_spos <= gc:
                break
    #end for gc in range(0,len(gene_coorDic)):
    return gene_spos,newpos
#def Get_gene_spos():


def Get_gene_spos_Reverse(g_pos,gs,rs,gene_coorDic,gene_cigarDic,gene_seqDic):
    diff = rs - gs
    basecigar = ''
    newcigar = ''
    newpos = 0
    gene_spos = 0
    seq_spos = 0
    
    coorflag = 0#0:number,1:'*'
    
    seq_num= 0
    for gc_c in range(0,len(gene_coorDic)):
        gc = len(gene_coorDic) - gc_c -1
        if gc == 0:
            gene_spos = gc
        else:
            gene_spos = gc+1
            
        
        if (gene_seqDic[gc] == '.') or (gene_seqDic[gc] == '*'):
            #print 'seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]',seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc],'----skip--------'
            continue
        else:
            #print 'seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]',seq_spos,gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc],'gene_spos:',gene_spos
            '''
            if gc == 0:
                gene_spos = gc
            else:
                gene_spos = gc-1
            '''    
                
                
                
            if seq_spos == diff:
                break
                    
            #end if seq_spos == diff:
            seq_num = gc
            seq_spos = seq_spos + 1
            #end if gene_spos == diff:
        # end else
        
    #end for gc in range(diff,len(gene_coorDic)):
    start_base = 0
    for gc in range(0,len(gene_coorDic)):
        if (gene_coorDic[gc]=='*') or (gene_cigarDic[gc]=='S'):
            continue
        else:
            newpos = g_pos+start_base
            start_base = start_base + 1
            if gene_spos == gc:
                break
    #end for gc in range(0,len(gene_coorDic)):
    return gene_spos,newpos
#def Get_gene_spos_Reverse():


def getNewCigar_Foward_Foward(g_pos,gs,rs,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,
                gene_seqDic,read_seqDic):
    diff = rs - gs
    basecigar = ''
    newcigar = ''
    newpos = 0
    seq_spos = 0
    
    coorflag = 0#0:number,1:'*'
   
    
    
    #print 'gene_spos:',gene_spos
    
    #print 'Cigar----------'
    for rc in range(0,len(read_coorDic)):
        #print rc,read_coorDic[rc],read_cigarDic[rc],read_seqDic[rc]
        if read_cigarDic[rc] == 'M':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'S':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == 'M':
        
        elif read_cigarDic[rc] == 'I':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'I':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'S'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
        #end if read_cigarDic[rc] == 'I':
        
        elif read_cigarDic[rc] == 'D':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
        #end if read_cigarDic[rc] == 'D':
        
        elif read_cigarDic[rc] == 'N':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
        #end if read_cigarDic[rc] == 'N':
        
        elif read_cigarDic[rc] == 'S':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'S'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
        #end if read_cigarDic[rc] == 'S':
        
        elif read_cigarDic[rc] == 'H':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
        #end if read_cigarDic[rc] == 'H':
        
        elif read_cigarDic[rc] == 'P':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
        #end if read_cigarDic[rc] == 'P':
        
        if read_cigarDic[rc] == '=':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == '=':
        
        if read_cigarDic[rc] == 'X':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == 'X':
        
    #end for rc in range(0,len(read_coorDic)):

    cigar_c=0#
    string_c = ''
    count = 0
    changeCigar = False
    for c in basecigar:
        count = count + 1
        if count == 1:
            string_c = c
        if string_c == c:
            changeCigar = False
            cigar_c = cigar_c + 1
        else:
            newcigar = newcigar + str(cigar_c) + string_c
            changeCigar = True
            string_c = c
            cigar_c = 1
        if count == len(basecigar):
            newcigar = newcigar + str(cigar_c) + string_c
    #end for c in basecigar:
    #print newcigar
    return newcigar
#end getNewCigar_Foward_Foward

def getNewCigar_Reverse_Foward(g_pos,gs,rs,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,
                gene_seqDic,read_seqDic):
    diff = rs - gs
    basecigar = ''
    newcigar = ''
    newpos = 0
    seq_spos = 0
    
    coorflag = 0#0:number,1:'*'
   
    
    
    #print 'gene_spos:',gene_spos
    
    #print 'Cigar----------'
    ss = 0
    for rc_c in range(0,len(read_coorDic)):
        
        rc = len(read_coorDic) - rc_c-1
        #print rc,read_coorDic[rc],read_cigarDic[rc],read_seqDic[rc]
        if read_cigarDic[rc] == 'M':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                
                
                
                if (rc == 0) and (gc == -1):
                    break
                
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                    ##print 'ss',ss
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'S':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == '=':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                
                
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == 'M':
        
        elif read_cigarDic[rc] == 'I':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'S'
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
        #end if read_cigarDic[rc] == 'I':
        
        elif read_cigarDic[rc] == 'D':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
        #end if read_cigarDic[rc] == 'D':
        
        elif read_cigarDic[rc] == 'N':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
        #end if read_cigarDic[rc] == 'N':
        
        elif read_cigarDic[rc] == 'S':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc - 1
                    ss = gc
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc - 1
                    ss = gc
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'S'
                    gene_spos = gc + 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc - 1
                    ss = gc
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
        #end if read_cigarDic[rc] == 'S':
        
        elif read_cigarDic[rc] == 'H':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
        #end if read_cigarDic[rc] == 'H':
        
        elif read_cigarDic[rc] == 'P':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    ss = gc
                    break
        #end if read_cigarDic[rc] == 'P':
        
        if read_cigarDic[rc] == '=':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    
                    #print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == '=':
        
        if read_cigarDic[rc] == 'X':
            for gc in range(gene_spos,len(gene_coorDic)):
                if rc == len(read_coorDic)-1:
                    ss = gc
                else:
                    
                    ##print 'gc,ss',gc,ss
                    if gc >= ss:
                        gc = ss-1
                        #print 'gc,ss',gc,ss
                        
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    ss = gc
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc - 1
                    ss = gc
                    break
                
                
                
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == 'X':
        
    #end for rc in range(0,len(read_coorDic)):

    cigar_c=0#
    string_c = ''
    count = 0
    changeCigar = False
    for c in basecigar:
        count = count + 1
        if count == 1:
            string_c = c
        if string_c == c:
            changeCigar = False
            cigar_c = cigar_c + 1
        else:
            newcigar = newcigar + str(cigar_c) + string_c
            changeCigar = True
            string_c = c
            cigar_c = 1
        if count == len(basecigar):
            newcigar = newcigar + str(cigar_c) + string_c
    #end for c in basecigar:
    #print newcigar
    return newcigar,newpos
#end getNewCigar_Reverse_Foward


def reverseSeq(g_seq):
    new_seq = ''
    for c in range(0,len(g_seq)):
        if g_seq[c:c+1] == 'A':
            new_seq = 'T'+new_seq
        elif g_seq[c:c+1] == 'T':
            new_seq = 'A'+new_seq
        elif g_seq[c:c+1] == 'G':
            new_seq = 'C'+new_seq
        elif g_seq[c:c+1] == 'C':
            new_seq = 'G'+new_seq
        else:
            new_seq = g_seq[c:c+1]+new_seq
        
    return new_seq
#end reverseSeq(g_seq)

def reverseCigar(g_cigar):
    new_cigar = ''
    num = ''
    
    for c in range(0,len(g_cigar)):
        numboolean = g_cigar[c:c+1].isdigit()
        if numboolean == True:
            num = num + g_cigar[c:c+1]
        else:
            new_cigar = num + g_cigar[c:c+1] + new_cigar
            num = ''
    return new_cigar
#end reverseCigar(g_cigar)

def reverse(r_seq):
    new_seq = ''
    for c in range(0,len(r_seq)):
        if r_seq[c:c+1] == 'A':
            new_seq = 'T'+new_seq
        elif r_seq[c:c+1] == 'T':
            new_seq = 'A'+new_seq
        elif r_seq[c:c+1] == 'G':
            new_seq = 'C'+new_seq
        elif r_seq[c:c+1] == 'C':
            new_seq = 'G'+new_seq
        else:
            new_seq = r_seq[c:c+1] + new_seq
    return new_seq
#end reverse(r_seq)

def ReverseReversegetNewCigar(reverse_seq,g_seq,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,
                gene_seqDic,read_seqDic):
    basecigar = ''
    newcigar = ''
    newpos = 0
    seq_spos = 0
    coorflag = 0#0:number,1:'*'
    for rc in range(0,len(read_coorDic)):
        if read_cigarDic[rc] == 'P':
            basecigar = basecigar + read_cigarDic[rc]
            continue
        elif read_cigarDic[rc] == 'I':
            basecigar = basecigar + read_cigarDic[rc]
            continue
        elif (read_cigarDic[rc] == 'D') and (read_cigarDic[rc] == 'N'):
            for gc in range(gene_spos,len(gene_coorDic)):
                if (gene_cigarDic[gc] == 'D') and (gene_cigarDic[gc] == 'N') and (gene_cigarDic[gc] == 'P'):
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    continue
                else:
                    if (gene_cigarDic[rc] == 'M') or (gene_cigarDic[rc] == '=') or (gene_cigarDic[rc] == 'X'):
                         basecigar = basecigar + read_cigarDic[rc]
                         gene_spos = gc + 1
                         break
                    else:
                         basecigar = basecigar + read_cigarDic[rc]
                         gene_spos = gc + 1
                         break
        else:
            #print rc,read_coorDic[rc],read_cigarDic[rc],read_seqDic[rc]
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'D':
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    continue
                elif gene_cigarDic[gc] == 'N':
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    continue
                elif gene_cigarDic[gc] == 'P':
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    continue
                else:
                    if (gene_cigarDic[rc] == 'M') or (gene_cigarDic[rc] == '=') or (gene_cigarDic[rc] == 'X'):
                         basecigar = basecigar + read_cigarDic[rc]
                         gene_spos = gc + 1
                         break
                    else:
                         basecigar = basecigar + gene_cigarDic[gc]
                         gene_spos = gc + 1
                         break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end else:
    #end for rc in range(0,len(read_coorDic)):
    
    cigar_c=0#
    string_c = ''
    count = 0
    changeCigar = False
    for c in basecigar:
        count = count + 1
        if count == 1:
            string_c = c
        if string_c == c:
            changeCigar = False
            cigar_c = cigar_c + 1
        else:
            newcigar = newcigar + str(cigar_c) + string_c
            changeCigar = True
            string_c = c
            cigar_c = 1
        if count == len(basecigar):
            newcigar = newcigar + str(cigar_c) + string_c
    #end for c in basecigar:
    return newcigar,newpos
#end def getNetCigar():


def ReversegetNewCigar(reverse_seq,g_seq,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,
                gene_seqDic,read_seqDic):
    
    basecigar = ''
    newcigar = ''
    newpos = 0
    seq_spos = 0
    
    coorflag = 0#0:number,1:'*'
    
    #print 'Cigar-----------------'
    for rc in range(0,len(read_coorDic)):
        #print rc,read_coorDic[rc],read_cigarDic[rc],read_seqDic[rc]
        if read_cigarDic[rc] == 'M':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'S':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == 'M':
        
        elif read_cigarDic[rc] == 'I':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'I':
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + 'S'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
        #end if read_cigarDic[rc] == 'I':
        
        elif read_cigarDic[rc] == 'D':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
        #end if read_cigarDic[rc] == 'D':
        
        elif read_cigarDic[rc] == 'N':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
        #end if read_cigarDic[rc] == 'N':
        
        elif read_cigarDic[rc] == 'S':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + 'S'
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + 'H'
                    gene_spos = gc + 1
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
        #end if read_cigarDic[rc] == 'S':
        
        elif read_cigarDic[rc] == 'H':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
        #end if read_cigarDic[rc] == 'H':
        
        elif read_cigarDic[rc] == 'P':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + read_cigarDic[rc]
                    break
        #end if read_cigarDic[rc] == 'P':
        
        if read_cigarDic[rc] == '=':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == '=':
        
        if read_cigarDic[rc] == 'X':
            for gc in range(gene_spos,len(gene_coorDic)):
                if gene_cigarDic[gc] == 'M':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'I':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'D':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'N':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    asecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'S':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'H':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == 'P':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    basecigar = basecigar + gene_cigarDic[gc]
                elif gene_cigarDic[gc] == '=':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
                elif gene_cigarDic[gc] == 'X':
                    #print '-------------',gc,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                    if gene_coorDic[gc] != '*':
                             newpos = int(gene_coorDic[gc])
                    basecigar = basecigar + gene_cigarDic[gc]
                    gene_spos = gc + 1
                    break
            #end for gc in range(gene_spos,len(gene_coorDic)):
        #end if read_cigarDic[rc] == 'X':
        
    #end for rc in range(0,len(read_coorDic)):
    
    
    cigar_c=0#
    string_c = ''
    count = 0
    changeCigar = False
    for c in basecigar:
        count = count + 1
        if count == 1:
            string_c = c
        if string_c == c:
            changeCigar = False
            cigar_c = cigar_c + 1
        else:
            newcigar = newcigar + str(cigar_c) + string_c
            changeCigar = True
            string_c = c
            cigar_c = 1
        if count == len(basecigar):
            newcigar = newcigar + str(cigar_c) + string_c
    #end for c in basecigar:
    #print '***************************************',newcigar
    return newcigar,newpos
#end def getNetCigar():



def MainSubCompare(new_read_file,rnameList,SAM_atline):
    if len(rnameList) == 1:
        f=open(new_read_file + "." + rnameList[0], "w")
    else:
        f=open(new_read_file, "w")

    for line in SAM_atline:
        f.write(line)
    #end for line in SAM_atline:
    
    line_ENTER = 0
    break_c = 0#break counter
    for ch in rnameList:
        '''
        if ch != 'scaffold_1048':
            continue
        '''
        
        log.info('Now converting:' + ch)
        rgdic = {} #dic[read_ID]=gene_ID
        gdic = {} #dic[gene_ID] = SAM_line
        gene_file = os.path.join(working_dir,'GeneOnMainRef.'+ch+'.sam')
        rgfile = os.path.join(working_dir, 'ReadOnGeneList.'+ch+'.tab')
        read_file = os.path.join(working_dir, 'MappedRead.'+ch+'.sam')
        #new_read_file = 'Output.'+ch+'.sam'
        geneSamList = [] #[sam_line]
        #rSamList = []
        id = ''
        write_line = 0
        
        gene_coorDic = {}
        gene_cigarDic = {}
        gene_seqDic = {}
        gene_flag = 0
        p_gid = ''
        p_gref = ''
        p_gpos = ''
        p_greverse = ''
        
        for line in open(rgfile):
            itemList = line[:-1].split('\t')
            rg_gID = itemList[0]
            gs = int(itemList[1])
            rg_gdirect = itemList[3]
            rg_rID = itemList[4]
            rs = int(itemList[5])
            read_sam  = itemList[6]
            
            '''
            if rg_rID != 'SRR072806.1534907':
                continue
            '''
            
            if rs == gs:
                continue
            
            #print line
            rID = itemList[7]
            r_cigar = itemList[12]
            r_seq = itemList[16]
                
            rreverse = reverseBoolean(itemList[8])#foward or reverse
            
            newcigar = ''#read new cigar
            newpos = 0
            
            if id != rg_gID:
                geneSamList = getSamIDList(ch,rg_gID, gene_file)
                id = rg_gID
                if len(geneSamList) == 0:
                    log.info('No mapping: (Read)' + rg_rID  + ',(Gene)' + itemList[0] )
                    continue
                
                
            for g in geneSamList:
                g_itemList = g[:-1].split('\t')
                gID = getsubString_a(g_itemList[0],'|')
                g_ref = g_itemList[2]
                g_pos = g_itemList[3]
                g_cigar = g_itemList[5]
                g_seq = g_itemList[9]
                
                if g_ref == '*':
                    continue
                
                
                greverse = reverseBoolean(g_itemList[1])#foward or reverse gs
                #log.info('Now converting: '+ rg_rID + " len:" + str(len(g_seq))+ " dir:" + str(rreverse) + ',' + itemList[0] + " dir:" + str(greverse))
                
                
                #print greverse
                if greverse == False:
                    if rg_gdirect == '+':
                        read_coorDic,read_cigarDic,read_seqDic = addCigarSeq(r_cigar,r_seq)
                        if not ((greverse == p_greverse) and 
                                (gID == p_gid) and (g_ref == p_gref) and (g_pos == p_gpos)):
                            p_gid = gID
                            p_gref = g_ref
                            p_gpos = g_pos
                            p_greverse = greverse
                            gene_coorDic = {}
                            gene_cigarDic = {}
                            gene_seqDic = {}
                            gene_coorDic,gene_cigarDic,gene_seqDic = addCigarSeq(g_cigar,g_seq)
                            sys.stderr.write("\n")
                            log.info('Now converting: '+ rg_gID + " len:" + str(len(g_seq))+ " dir:" + str(rreverse) + "(" + gID + "," + g_ref+ ","+ str(g_pos) + "," + str(greverse)+ ")")
                        else:
                            sys.stderr.write(".")

                        #end if (gID != p_gid) and (g_ref != p_gref) and (g_pos != p_gpos):
                        gene_spos,newpos = Get_gene_spos(int(g_pos),gs,rs,gene_coorDic,gene_cigarDic,gene_seqDic)
                        #print 'gene_spos,newpos',gene_spos,newpos
                        newcigar = getNewCigar_Foward_Foward(int(g_pos),gs,rs,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,gene_seqDic,read_seqDic)
                    #end if rg_gdirect == '+':
                    
                    elif rg_gdirect == '-':
                        
                        reverse_rseq = reverseSeq(r_seq)
                        #print 'reverse_rseq:', reverse_rseq
                        reverse_rcigar = reverseCigar(r_cigar)
                        read_coorDic,read_cigarDic,read_seqDic = addCigarSeq(reverse_rcigar,reverse_rseq)
                        if not ((greverse == p_greverse) and 
                                (gID == p_gid) and (g_ref == p_gref) and (g_pos == p_gpos)):
                            p_gid = gID
                            p_gref = g_ref
                            p_gpos = g_pos
                            p_greverse = greverse
                            gene_coorDic = {}
                            gene_cigarDic = {}
                            gene_seqDic = {}
                            gene_coorDic,gene_cigarDic,gene_seqDic = addCigarSeq(g_cigar,g_seq)
                            sys.stderr.write("\n")
                            log.info('Now converting: '+ rg_gID + " len:" + str(len(g_seq))+ " dir:" + str(rreverse) + "(" + gID + "," + g_ref+ ","+ str(g_pos)+ "," + str(greverse)+ ")")
                        else:
                            sys.stderr.write(".")

                        #end if (gID != p_gid) and (g_ref != p_gref) and (g_pos != p_gpos):
                        
                        gene_spos,newpos = Get_gene_spos_Reverse(int(g_pos),gs,rs,gene_coorDic,gene_cigarDic,gene_seqDic)
                        newcigar, newpos = getNewCigar_Reverse_Foward(int(g_pos),gs,rs,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,gene_seqDic,read_seqDic)
                        #print 'newpos:',newpos
                        
                        newcigar = reverseCigar(newcigar)
                        str_newpos = str(newpos)
                        position = int(g_pos)
                        seq_num = 0
                        change_pos = 0
                        
                        
                        for gc in range(0,len(gene_coorDic)):
                            if (gene_cigarDic[gc] == 'S') or (gene_cigarDic[gc] == 'I') or (gene_cigarDic[gc] == 'P') or (gene_cigarDic[gc] == 'H'):
                                if (gene_cigarDic[gc] == 'S') or (gene_cigarDic[gc] == 'I'):
                                    #print 'position,seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]','*',seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                                    if str_newpos == gene_coorDic[gc]:
                                        
                                        for gc_c in range(gc,len(gene_coorDic)):
                                            ##print position,seq_num,gene_coorDic[gc_c],gene_cigarDic[gc_c],gene_seqDic[gc_c]
                                            if (gene_cigarDic[gc_c] == 'S') or (gene_cigarDic[gc_c] == 'I') or (gene_cigarDic[gc_c] == 'P') or (gene_cigarDic[gc_c] == 'H'):
                                                continue
                                            else:
                                                change_pos = position
                                                #print 'chage_pos:', change_pos
                                                break
                                        #end for gc_c in range(gc,len(gene_coorDic)):
                                    #if str_newpos == gene_coorDic[gc]:
                                    seq_num = seq_num + 1
                                
                                if str_newpos == gene_coorDic[gc]:
                                    change_pos = position
                                    break
                                continue
                            else:
                                #print 'position,seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]',position,seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                                
                                if  (gene_cigarDic[gc] == 'M') or (gene_cigarDic[gc] == '=') or (gene_cigarDic[gc] == 'X'):
                                    seq_num = seq_num + 1
                                if str_newpos == gene_coorDic[gc]:
                                    change_pos = position
                                    break
                                position = position + 1
                    #end for gc_c in range(0,len(gene_coorDic)):
                        
                        
                    #end elif rg_gdirect == '-':
                    
                    if line_ENTER  != 0:
                        f.write('\n')
                    #end(if line_ENTER  != 0:)
                    line_ENTER  = line_ENTER + 1
                    
                    tab_c = 0
                    for tab in range(7,len(itemList)):
                        if tab_c == 2:
                            f.write(g_ref)
                        elif tab_c == 3:
                            if rg_gdirect == '+':
                                f.write(str(newpos))
                            elif rg_gdirect == '-':
                                f.write(str(change_pos))
                        elif tab_c == 5:
                            f.write(newcigar)
                        elif tab_c == 9:
                            if rg_gdirect == '+':
                                f.write(r_seq)
                            else:
                                f.write(reverse_rseq)
                        else:
                            f.write(itemList[tab])
                        tab_c = tab_c + 1
                        if tab_c != len(itemList)-7:
                            f.write("\t")
                    #end for tab in r_itemList:
                #end(if greverse == False:)
                
                
                if greverse == True:
                    if rg_gdirect == '+':
                        read_coorDic,read_cigarDic,read_seqDic = addCigarSeq(r_cigar,r_seq)
                        if not ((greverse == p_greverse) and 
                                (gID == p_gid) and (g_ref == p_gref) and (g_pos == p_gpos)):
                            reverse_seq = reverseSeq(g_seq)
                            reverse_cigar = reverseCigar(g_cigar)
                            p_gid = gID
                            p_gref = g_ref
                            p_gpos = g_pos
                            p_greverse = greverse
                            gene_coorDic = {}
                            gene_cigarDic = {}
                            gene_seqDic = {}
                            gene_coorDic,gene_cigarDic,gene_seqDic = addCigarSeq(reverse_cigar,reverse_seq)
                            sys.stderr.write("\n")
                            log.info('Now converting: '+ rg_gID + " len:" + str(len(g_seq))+ " dir:" + str(rreverse) + "(" + gID + "," + g_ref+ ","+ str(g_pos)+ "," + str(greverse)+ ")")
                        else:
                            sys.stderr.write(".")
                        #end if (gID != p_gid) and (g_ref != p_gref) and (g_pos != p_gpos):
                        
                        
                        gene_spos,newpos = Get_gene_spos(int(g_pos),gs,rs,gene_coorDic,gene_cigarDic,gene_seqDic)
                        newcigar,newpos = ReversegetNewCigar(reverse_seq,g_seq,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,gene_seqDic,read_seqDic)
                        
                        #print 'newpos:', newpos
                        newcigar = reverseCigar(newcigar)
                        newr_seq =reverse(r_seq)
                        str_newpos = str(newpos)
                        
                        position = int(g_pos)
                        seq_num = 0
                        change_pos = 0
                        
                        for gc_c in range(0,len(gene_coorDic)):
                            gc = len(gene_coorDic) - gc_c -1
                            
                            if (gene_cigarDic[gc] == 'S') or (gene_cigarDic[gc] == 'I') or (gene_cigarDic[gc] == 'P') or (gene_cigarDic[gc] == 'H'):
                                if (gene_cigarDic[gc] == 'S') or (gene_cigarDic[gc] == 'I'):
                                    #print 'position,seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]','*',seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                                    if str_newpos == gene_coorDic[gc]:
                                        for gc_c in range(gc,len(gene_coorDic)):
                                            if (gene_cigarDic[gc_c] == 'S') or (gene_cigarDic[gc_c] == 'I') or (gene_cigarDic[gc_c] == 'P') or (gene_cigarDic[gc_c] == 'H'):
                                                continue
                                            else:
                                                change_pos = position
                                                #print 'chage_pos:', change_pos
                                                break
                                        ##end for gc_c in range(gc,len(gene_coorDic)):
                                    #if str_newpos == gene_coorDic[gc]:
                                    seq_num = seq_num + 1
                                if str_newpos == gene_coorDic[gc]:
                                    change_pos = position
                                    break
                                continue
                            #end if (gene_cigarDic[gc] == 'S') or (gene_cigarDic[gc] == 'I') or (gene_cigarDic[gc] == 'P') or (gene_cigarDic[gc] == 'H'):
                            else:
                                #print 'position,seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]',position,seq_num,gene_coorDic[gc],gene_cigarDic[gc],gene_seqDic[gc]
                                
                                if  (gene_cigarDic[gc] == 'M') or (gene_cigarDic[gc] == '=') or (gene_cigarDic[gc] == 'X'):
                                    seq_num = seq_num + 1
                                if str_newpos == gene_coorDic[gc]:
                                    change_pos = position
                                    break
                                position = position + 1
                        #end for gc_c in range(0,len(gene_coorDic)):
                        
                        
                        if line_ENTER  != 0:
                            f.write('\n')
                        line_ENTER  = line_ENTER + 1
                        
                        tab_c = 0 
                        for tab in range(7,len(itemList)):
                            if tab_c == 2:
                                f.write(g_ref)
                            elif tab_c == 1:
                                int_flag = int(itemList[8])
                                if rreverse == False:
                                    f.write(str(int_flag | 0b10000))
                                elif rreverse == True:
                                    f.write(str(int_flag ^ 0b10000))
                            elif tab_c == 3:
                                f.write(str(change_pos))
                            elif tab_c == 5:
                                f.write(newcigar)
                            elif tab_c == 9:
                                f.write(newr_seq)
                            else:
                                f.write(itemList[tab]) 
                            tab_c = tab_c + 1
                            if tab_c != len(itemList)-7:
                                f.write("\t")
                        #end for tab in range(7,len(itemList)):
                        
                        
                    
                    #end if rg_gdirect == '+':
                    
                    elif rg_gdirect == '-':
                        read_coorDic,read_cigarDic,read_seqDic = addCigarSeq(r_cigar,r_seq)
                        if not ((gID == p_gid) and (g_ref == p_gref) and (g_pos == p_gpos)):
                            p_gid = gID
                            p_gref = g_ref
                            p_gpos = g_pos
                            gene_coorDic = {}
                            gene_cigarDic = {}
                            gene_seqDic = {}
                            gene_coorDic,gene_cigarDic,gene_seqDic = addCigarSeq(g_cigar,g_seq)
                            sys.stderr.write("\n")
                            log.info('Now converting: '+ rg_gID + " len:" + str(len(g_seq))+ " dir:" + str(rreverse) + "(" + gID + "," + g_ref+ ","+ str(g_pos)+ "," + str(greverse)+ ")")
                        else:
                            sys.stderr.write(".")
                        #end if (gID != p_gid) and (g_ref != p_gref) and (g_pos != p_gpos):
                        
                        gene_spos,newpos = Get_gene_spos(int(g_pos),gs,rs,gene_coorDic,gene_cigarDic,gene_seqDic)
                        newcigar = getNewCigar_Foward_Foward(int(g_pos),gs,rs,gene_spos,gene_coorDic,gene_cigarDic,read_coorDic,read_cigarDic,gene_seqDic,read_seqDic)
                   
                        #print 'newpos:', newpos
                        
                        if line_ENTER  != 0:
                            f.write('\n')
                        #end(if line_ENTER  != 0:)
                        line_ENTER  = line_ENTER + 1
                        
                        
                        tab_c = 0
                        for tab in range(7,len(itemList)):
                            if tab_c == 2:
                                f.write(g_ref)
                            elif tab_c == 3:
                                f.write(str(newpos))
                            elif tab_c == 5:
                                f.write(newcigar)
                            elif tab_c == 9:
                                f.write(r_seq) 
                            else:
                                f.write(itemList[tab])
                            tab_c = tab_c + 1
                            if tab_c != len(itemList)-7:
                                f.write("\t")
                        #end for tab in range(7,len(itemList)):
                    #end elif rg_gdirect == '-':
                
                #end(if (greverse == True) and (rreverse == False):)
                
                else:
                    continue
                
                
                
                
             
    f.close()
    #print 'Make file:',new_read_file
    
#end (MainSubCompare(out_file,rnameList,SAM_atline))



def getGeneList(gene_sam_file):
    id_dic = {} # {id,number list}
    sam_line_dic = {} #{number, sam_line}
    
    log.info('Get GeneList')
    fileline = 0
    try:#file open
        for line in open(gene_sam_file):
            itemList = line[:-1].split('\t')
            if itemList[0][:1] == '@':
                continue
            id = itemList[0][:itemList[0].find('|')]
            
            if (id in id_dic) == False:
                #print id,fileline
                id_dic [id] = str(fileline)
            #end if (id in id_dic) == False:
                
            elif (id in id_dic) == True:
                value = id_dic[id]
                value = value + ',' + str(fileline)
                id_dic[id] = value
            #end elif (id in id_dic) == True:
            
            
            sam_line_dic[str(fileline)] = line
            fileline = fileline + 1
            '''
            if fileline == 5000:
                break
            '''
        #end for line in open(rgfile):
    #end try:#file open
    except IOError, (errno, msg):
        print 'except: Cannot open "%s"' % gene_sam_file
        print '  errno: [%d] msg: [%s]' % (errno, msg)
    #end except IOError, (errno, msg):
    
    '''
    print 'id_dic:',len(id_dic)
    print 'sam_line_dic:',len(sam_line_dic)
    '''
    
    return id_dic,sam_line_dic
#end getGeneList(gene_sam_file)

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
    parser.add_option("-n", "--number", action="store", dest="number",
                      default="-1")

    (opt, args) = parser.parse_args()

    number = int(opt.number)

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

    config_file = opt.config
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    global working_dir
    try:
        log.info("Importing settings from %s" % config_file )
        read_sam_file = config.get("combine_reads", "read_sam_file")
        log.info("read_sam_file=" + read_sam_file)
        gene_sam_file = config.get("removemultiple", "gene_seq_uniq_sam")
        log.info("gene_sam_file=" + gene_sam_file)
        out_file = config.get("output", "final_sam_file")
        log.info("out_file" + out_file)
        working_dir = config.get("combine_reads", "working_dir")
        log.info("working_dir=" + working_dir)

        # To show next command
        original_fasta_file = config.get("global", "original_fasta_file")
        target_fasta_file = config.get("global", "target_fasta_file")
        read_sam_file=config.get("combine_reads", "read_sam_file")
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

    rnameList = []
    SAM_header = []
    rnameList = getRefName(read_sam_file,rnameList)
    getChrSAMfile(gene_sam_file,rnameList)
    SAM_header = getSAMheader(gene_sam_file)
    if number == -1:
        MainSubCompare(out_file,rnameList,SAM_header)
    else:
        MainSubCompare(out_file,[rnameList[number]],SAM_header)

    p = re.compile("(.+).[sS][aA][mM]$")
    m = p.match(out_file)
    out_bam = out_file
    if m:
        out_bam = m.group(1)

    print("\n=====")
    print("Finished the coordinate convert. Final SAM file is %s." % out_file)
    print("To visualize the file in IGV, you have to convert the SAM file into BAM file.")
    print("Eg.")
    print("$ samtools view -bt %s -o your_read.bam %s" % (target_fasta_file, out_file))
    print("$ samtools sort your_read.bam %s" % out_bam)
    print("$ samtools index %s.bam" % out_bam)
    print("After this, you import %s.bam (and %s.bai) into IGV." % (out_bam,out_bam))
    
#end(def main():)



if __name__ == "__main__": 
    main()
