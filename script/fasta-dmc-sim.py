#! /usr/bin/env python
# 
# Read FASTA-format DNA sequences with symbols indicating
# methylation rates, and randomly assign differentially
# methylated cytosines.
# 
# Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/

import fileinput, optparse, os, random, signal, sys

def fastaInput(lines):
    title = ""
    for line in lines:
        line = line.rstrip()
        if line.startswith(">"):
            if title: yield title, seq
            seq = []
            title = line
        else:
            seq.extend(line)
    if title: yield title, seq

def writeFasta(title, seq):
    linesize = 50
    seqlen = len(seq)
    beg = 0
    print title
    while beg < seqlen:
        end = beg + linesize
        print "".join(seq[beg:end])
        beg = end

def dmcSeq(seq, dmcIndex):
    for dmc in dmcIndex:
        oldBase = seq[dmc[0]]
        if oldBase in "tvdc": # forward strand
            otherBases = [i for i in "tvdcC" if i != oldBase]
        else:                 # reverse strand
            otherBases = [i for i in "abhgG" if i != oldBase]
        newBase = random.choice(otherBases)
        seq[dmc[0]] = newBase
        dmc[0] = str(dmc[0])
        dmc[2] += ","+newBase

def fastaDmcSim(opts, args):
    dmcRate = float(args[0])
    random.seed(1414)  # seed
    if opts.f: f = open(opts.f, "w")
    for title, seq in fastaInput(fileinput.input(args[1])):
        seqlen = len(seq)
        methylIndex = [[i, str(i+1), seq[i]] for i in range(seqlen) if seq[i] in "tvdcabhg"]

        dmcNum = int(len(methylIndex) * dmcRate)
        dmcIndex = random.sample(methylIndex, dmcNum)
        dmcIndex.sort()
        seq, dmcIndex = dmcSeq(seq, dmcIndex)
        writeFasta(title, seq)
        if os.path.exists(opts.f):
            for dmc in dmcIndex:
                f.write("\t".join([title.lstrip(">")]+dmc) + "\n")
    if os.path.exists(opts.f):
        f.close()


if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    usage = "%prog dmcFraction fastaFile(s)"
    description = "Randomly assign differentially methylated cytosines (DMC) with specific rates."

    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-f", metavar="FILE", help="write DMC positions to bed FILE")
    opts, args = op.parse_args()

    if len(args) != 2: op.error("I need 2 file names")
    
    try: fastaDmcSim(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
