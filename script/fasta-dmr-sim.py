#!/usr/bin/env python
#
# Read FASTA-format DNA sequences with symbols indicating
# methylation rates, and randomly assign differentially
# methylated cytosines.
# fasta-dmr-sim
#
# Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/

import fileinput
import optparse
import os
import random
import signal
import sys
import re

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
    return

def writeFasta(title, seq):
    linesize = 50
    seqlen = len(seq)
    beg = 0
    print title
    while beg < seqlen:
        end = beg + linesize
        print "".join(seq[beg:end])
        beg = end
    return

def dmrSeq(seqName, seq, context, dmcIndex, numDmr, lenDmr):
    # len(dmcIndex) is far more than numDmr * lenDmr
    frac = numDmr * lenDmr / len(dmcIndex)
    if frac > 0.9:
        raise "DMR parameters (numDmr and lenDmr) are too large."
        return
    dmrs = []
    size = len(dmcIndex)
    index = range(size)
    while len(dmrs)!=numDmr:
        i = random.choice(index)
        if len(dmrs)==0 or (i!=-1 and (i+lenDmr)!=-1):
            index[i] = -1
            j = i + 1
            while dmcIndex[j] < dmcIndex[i]+lenDmr and j < size:
                index[j] = -1
                j += 1
            upDown = 'UP'
            if random.random() < 0.5: upDown = 'DOWN'
            dmrs[seqName, dmcIndex[i], dmcIndex[j], upDown, context]
    dmrs.sort()
    return dmrs

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
    return (seq, dmc)

def makeMethylIndex(seq, context):
    seqlen = len(seq)
    if context == 'ALL':
        methylIndex = [[i, str(i+1), seq[i]] for i in range(seqlen) if seq[i] in "tvdcabhg"]
    return methylIndex

def fastaDmcSim(opts, args):
    dmcRate = float(args[0])
    random.seed(1414)  # seed
    if opts.f:
        f = open(opts.f, "w")
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
    return

def normalizeContext(opt_c="ALL"):
    tmp = '|'.split(opt_c)
    flag = 0
    for i in tmp:
        if re.match('(cg|cpg)', i, re.IGNORECASE):
            flag = flag | 0x1
        elif re.match('(chg|cphpg)', i, re.IGNORECASE):
            flag = flag | 0x2
        elif re.match('(chh|cphph)', i, re.IGNORECASE):
            flag = flag | 0x4
        else:
            print >>sys.stderr, "Warn: \"%s\" is not a valid context keyword." % i
    set_str = ['', '', '']
    context = 'ALL'
    if flag | 0x1:
        set_str[0] = 'CpG'
    elif flag | 0x2:
        set_str[1] = 'CHG'
    elif flag | 0x4:
        set_str[2] = 'CHH'
    if flag != 0x7:
        context = '|'.join(set_str)
    return context

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    usage = "%prog [OPTIONS] <dmcFraction> <fastaFile> [<fastaFile> ...]"
    description = "Randomly assign differentially methylated cytosines (DMC) with specific rates."

    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-f", metavar="FILE", help="Write DMC positions to bed FILE")
    op.add_option("-c", metaval="STRING", help="Specify DMC context (eg. CG, CG|CHG, CG|CHH, ALL) default ALL")
    opts, args = op.parse_args()

    if len(args) != 2:
        op.error("I need 2 file names")

    context = normalizeContext(opt.c)

    try:
        fastaDmcSim(opts, args)
    except KeyboardInterrupt:
        pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
    sys.exit(0)
