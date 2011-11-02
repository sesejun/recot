#!/usr/bin/env python
"""Extract genes and their reguratory region sequences.
Configurations have to be described in ``settings.ini''

"""
__author__ =  'Akiko Izawa and Jun Sese'
__version__=  '0.1'

import ConfigParser
import logging
import sys, os
import re
from optparse import OptionParser
def cigerLength(ciger):
    #print 'ciger:',ciger
    num=''
    ciger_len = 0
    
    
    for c in ciger:
        if c == '*':
            ciger_len = 0
            break
        if num == '':
            num = num + c
        elif c.isdigit() == True:
            num = num + c
        elif c.isdigit() == False:
            #print 'c=', c, ',num=',num
            if c == 'M' or c == 'I' or c == 'S' or c == '=' or c == 'X':
                #print '--------ciger_len:', ciger_len ,'+',int(num), '=', ciger_len + int(num)
                ciger_len = ciger_len + int(num)
                
            num = ''
        
        
    #print '----->  ciger_len:',ciger_len
    return ciger_len
#end(def cigerLength(ciger))

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
            ciger = itemList[5]
            seq = itemList[9]
            
            ciger_length = cigerLength(ciger)
            
            if ciger == '*':
                f_out.write(line)
            elif len(seq) != ciger_length:
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
    #print 'Check ciger and seq'
    
    #print 'file:',home+file_name
    log.info("Input File: " + args[0])
    log.info("Output File: " + args[1])
    log.info("Error File: " + args[2])
    inputfile = args[0]
    outfile = args[1]
    errfile = args[2]
    #file='Pan_troglodytes.CHIMP2.1.64.dna.chromosome.all.geneSeq.on.Homo.all.sam'
    #outfile=home+'check_cigerandseq.sam'
    #error_file=home+'error.sam'
    loadFile(inputfile,outfile,errfile)
    print 'script end'
    

#end(def main():)


if __name__ == "__main__": 
    main()
