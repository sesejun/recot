#!/usr/bin/env python

"""join result sam files
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

import re, sys, errno, os
import glob
import logging
import ConfigParser
from optparse import OptionParser


def samJoin(working_dir, samname):
    fullpath = samname + ".*"
    log.info(fullpath)
    files = glob.glob(fullpath)
    p = re.compile(samname + "\.(.+)$")
    seqlist = []
    samfp = open(samname, "w")

    log.info("Reading headers")
    for file in files:
        for line in open(file, "r"):
            if line[0] == "@":
                seqlist.append(line)
            else:
                break
    seqlist = list(set(seqlist))
    for line in sorted(seqlist):
        samfp.write(line)

    log.info("Concatinating SAM entries")
    for file in files:
        log.debug(file)
        for line in open(file, "r"):
            if line[0] != "@":
                if line[-1] == '\n':
                    samfp.write(line)
                else:
                    samfp.write(line + "\n")
    
    #m = p.match(out_file)
    #out_bam = out_file
    #if m:
    #    out_bam = m.group(1)


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

    config_file = opt.config
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    global working_dir
    try:
        out_file = config.get("output", "final_sam_file")
        log.info("out_file" + out_file)
        working_dir = config.get("combine_reads", "working_dir")
        log.info("working_dir=" + working_dir)

        # to show next command
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

    samJoin(working_dir, out_file)

    p = re.compile("(.+).[sS][aA][mM]$")
    m = p.match(out_file)
    out_bam = out_file
    if m:
        out_bam = m.group(1)

    print("\n=====")
    print("Finished concatination of sam files. Final SAM file is %s." % out_file)
    print("To visualize the file in IGV, you have to convert the SAM file into BAM file.")
    print("Eg.")
    print("$ samtools view -bt %s -o your_read.bam %s" % (target_fasta_file, out_file))
    print("$ samtools sort your_read.bam %s" % out_bam)
    print("$ samtools index %s.bam" % out_bam)
    print("After this, you import %s.bam (and %s.bai) into IGV." % (out_bam,out_bam))


if __name__ == "__main__": 
    main()
