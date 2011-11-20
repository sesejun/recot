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
def cigarLength(cigar):
    #print 'cigar:',cigar
    num=''
    cigar_len = 0
    
    
    for c in cigar:
        if c == '*':
            cigar_len = 0
            break
        if num == '':
            num = num + c
        elif c.isdigit() == True:
            num = num + c
        elif c.isdigit() == False:
            #print 'c=', c, ',num=',num
            if c == 'M' or c == 'I' or c == 'S' or c == '=' or c == 'X':
                #print '--------cigar_len:', cigar_len ,'+',int(num), '=', cigar_len + int(num)
                cigar_len = cigar_len + int(num)
                
            num = ''
        
        
    #print '----->  cigar_len:',cigar_len
    return cigar_len
#end(def cigarLength(cigar))

def loadFile(file,outfile,error_file):
    
    f_out = open(outfile, "w")
    f_error = open(error_file, "w")
    fileline = 0
    for line in open(file):
        itemList = line[:-1].split('\t')
        fileline = fileline + 1
        
        if len(itemList) < 11:
            if itemList[0][0:1] != '@':
                log.debug(line)
                f_error.write(line)
                continue
            # end if itemList[0][0:1] != '@':
            else:
                f_out.write(line)
                continue
            # end else:
        # end if len(itemList) < 11:
        else:
            id = itemList[0]
            flag = int(itemList[1])
            cigar = itemList[5]
            seq = itemList[9]
            
            cigar_length = cigarLength(cigar)
            
            if cigar == '*':
                f_out.write(line)
            elif len(seq) != cigar_length:
                log.debug(line)
                f_error.write(line)
            else:
                f_out.write(line)
        # end else:
            
        
        
    #end (for line in open(gff_file):)
    f_out.close()
    f_error.close()
#end def loadFile(file)

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
    logging.basicConfig(level=logging.DEBUG,
                        format="%(name)s: %(levelname)s @ %(asctime)s %(message)s")

    #home='/home/akiko/TestScript/'
    #file_name='A.lyrata.testfile.on.TAIR10.sam'
    #print 'Check cigar and seq'
    
    #print 'file:',home+file_name
    log.info("Input File: " + args[0])
    log.info("Output File: " + args[1])
    log.info("Error File: " + args[2])
    inputfile = args[0]
    outfile = args[1]
    errfile = args[2]
    #file='Pan_troglodytes.CHIMP2.1.64.dna.chromosome.all.geneSeq.on.Homo.all.sam'
    #outfile=home+'check_cigarandseq.sam'
    #error_file=home+'error.sam'
    loadFile(inputfile,outfile,errfile)
    print 'script end'
    

#end(def main():)


if __name__ == "__main__": 
    main()
