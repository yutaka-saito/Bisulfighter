#!/usr/bin/env python
"""
Bisulfighter::bsf-call

Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
by National Institute of Advanced Industrial Science and Technology (AIST)
is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
http://creativecommons.org/licenses/by-nc-sa/3.0/
"""

__version__= "0.9"

import sys
import os
import glob
import threading
import subprocess
import Queue
from datetime import datetime
import hashlib
from string import maketrans
import gzip
import logging
from time import sleep
from shutil import copy
import re
import pysam

class BsfCallBase(object):
    """
    base class for BsfCall, LastExecutor, McDetector.
    define functions that is used in each sub classes.
    """

    def splitFilePath(self, filePath):
        """
        split file path to directory path, file name (with file extension),
        file name (without file extension) and file extension.
        if extension of filePath is '.gz', '.gz' extension is ignored.
        """

        if self.isGzippedFile(filePath):
            dir_name, file_name = os.path.split(filePath[0:-3])
        else:
            dir_name, file_name = os.path.split(filePath)
        base_name, ext = os.path.splitext(file_name)
        if len(ext) > 1:
            ext = ext[1:]

        return (dir_name, file_name, base_name, ext)


    def readNameByReadFile(self, readFilePath):
        """
        get read name by read file path.
        """

        dir_name, file_name, read_name, ext = self.splitFilePath(readFilePath)

        return read_name


    def secondReadFilePathByFirstReadFilePath(self, readFile, secondReadType = None):
        """
        get second read file path by first read file path.
        if first read file path is '/path/to/read_file_1.fastq' second read file
        path is '/path/to/read_file_2.fastq'
        if secondReadType is specified, the extension of second read file is its
        value.
        """

        fpath = ""

        if self.isGzippedFile(readFile):
            dir_name, file_name = os.path.split(readFile[0:-3])
            basename, ext = os.path.splitext(file_name)
            if secondReadType:
                ext = ".%s" % secondReadType
            fpath = "%s/%s2%s.gz" % (dir_name, basename[0:-1], ext)
        else:
            dir_name, file_name = os.path.split(readFile)
            basename, ext = os.path.splitext(file_name)
            if secondReadType:
                ext = ".%s" % secondReadType
            fpath = "%s/%s2%s" % (dir_name, basename[0:-1], ext)

        return fpath


    def pairedEndReadNumbers(self):
        return (1, 2)


    def clearGap(self, seq):
        return seq.replace("-", "")


    def complementStartPosition(self, genomeLen, subseqStart, subseqLen):
        return genomeLen - subseqStart - subseqLen


    def complementSeq(self, seq):
        return seq.translate(maketrans("ATGCatgc", "TACGtacg"))[::-1]


    def mcContextType(self, genomeSeq, cBasePos):
        """
        get mC context type (CG, CHG, CHH) by genome sequence and C base position.
        if no mC context found, return None.
        """

        try:
            if genomeSeq[cBasePos + 1] == "G":
                return "CG"
            else:
                if genomeSeq[cBasePos + 2] == "G":
                    return "CHG"
                else:
                    return "CHH"
            return None

        except IndexError:
            return None


    def chrSort(self, a, b):
        return cmp(a, b)


    def strands(self):
        return ("+", "-")


    def mcContextTypes(self):
        return ("CG", "CHG", "CHH")


    def gzipFile(self, filePath, wait = True, log = False):
        """
        gzip file. If wait argument is False, without waiting for gzip process to be
        completed, this function returns immediately.
        """

        if log:
            logging.info("gzip start: %s" % filePath)

        dirpath, fname = os.path.split(filePath)
        cmd = "gzip %s" % fname
        p = subprocess.Popen(cmd, shell = True, cwd = dirpath)
        if wait:
            p.wait()

        if log:
            logging.info("gzip done: %s" % filePath)


    def isGzippedFile(self, filePath):
        return filePath[-3:] == ".gz"


    def isMafFile(self, filePath):
        f = open(filePath, "rb")
        data = f.read(1)
        f.close()
        if re.match("[\x20-\x7E]", data):
            f = open(filePath, "r")
            first_line = f.readline()
            if first_line[0:5] != "track" and first_line[0] != "#" and first_line[0:2] != "a ":
                f.close()
                return False

            if first_line[0:2] == "a ":
                cnt = 2
            else:
                cnt = 1

            while True:
                line = f.readline()
                if line == "":
                    break
                if line[0] == "#" or line.strip() == "":
                    continue

                if cnt == 1 and line[0:2] != "a ":
                    f.close()
                    return False

                if cnt == 2 and line[0:2] != "s ":
                    f.close()
                    return False

                if cnt == 3:
                    f.close()
                    if line[0:2] == "s ":
                        return True
                    else:
                        return False

                cnt += 1

            f.close()
            return False
        else:
            return False


    def isBamFile(self, filePath):
        bgzf_magic = b"\x1f\x8b\x08\x04"

        f = open(filePath, "rb")
        data = f.read(4)
        f.close()

        return data == bgzf_magic


    def isSamFile(self, filePath):
        f = open(filePath, "rb")
        data = f.read(1)
        if data == "@":
            tag = f.read(2)
            if tag == "HD" or tag == "SQ" or tag == "RG" or tag == "CO":
                f.close()
                return True

        f.seek(0)
        data = f.read(1)
        f.close()
        if re.match("[\x20-\x7E]", data):
            f = open(filePath, "r")
            line = f.readline()
            f.close()
            return len(line.split("\t")) > 10
        else:
            return False


    def scriptDir(self):
        return os.path.dirname(os.path.abspath(sys.argv[0]))


    def chrnoFromFastaDescription(self, description):
        """
        get chromosome number by fasta description line.
        """

        return description.strip()[1:].strip()


    def chrsByRefGenome(self, refGenome):
        """
        get all chromosome numbers that is included in the reference genome.
        """

        chrs = []
        for line in open(refGenome, "r"):
            line = line.strip()
            if line[0] == ">":
                chrs.append(self.chrnoFromFastaDescription(line))

        return chrs


    def oneChrGenomeSeq(self, refGenome, chrNo):
        """
        get genome sequence of specified chrmosome number.
        """

        in_target_chr = False
        buf = []
        for line in open(refGenome, 'r'):
            line = line.strip()
            if line[0] == ">":
                if in_target_chr:
                    break
                else:
                    cur_chr = self.chrnoFromFastaDescription(line)
                    if cur_chr == chrNo:
                        in_target_chr = True
            else:
                if in_target_chr:
                    buf.append(line)

        return "".join(buf)


    def lastalOpts(self, lastOpt):
        """
        get options for lastal by bsf-call --last option
        """

        return " ".join(lastOpt.split(","))


    def mergeOpts(self):
        return ""


    def filterOpts(self, mismapProb, scoreThres, isPairedEnd):
        """
        get filtering option. this option is specified to last-map-probs or
        last-pair-probs.
        """

        option = ""
        if isPairedEnd:
            option = "-m%s" % str(mismapProb)
        else:
            option = "-s%d -m%s" % (scoreThres, str(mismapProb))

        return option


    def isPairedEnd(self, readAttr):
        return self.pairedEndReadNumbers()[1] in readAttr


    def jobIdByQsubOutput(self, qsubOutput):
        """
        get job id submitted by qsub command and its output.
        """

        fields = qsubOutput.strip().split()

        return fields[2]


    def waitForSubmitJobs(self, jobIds, checkInterval = 10):
        """
        wait for all jobs that have been submitted with qsub command to finish.
        """

        error_jobs = []
        while True:
            all_done = True
            qstat = os.popen("qstat")
            for line in qstat:
                fields = line.strip().split()
                if fields[0] in jobIds:
                    if fields[4].find('E') > 0:
                        if fields[0] not in error_jobs:
                            logging.fatal("Error has occurred: Job ID=%s" % fields[0])
                            error_jobs.append(fields[0])
                    else:
                        all_done = False
                        break
            qstat.close()

            if all_done:
                break
            else:
                sleep(checkInterval)

        return


    def logJobSubmit(self, msg, jobId, cmd = None):
        logging.info("Submit job: %s --> Job ID = %s" % (msg, jobId))
        if cmd:
            logging.info(cmd)


    def bamMapq2Mismap(self, mapq):
        return pow(0.1, (float(mapq) / 10))


    def getAllMappingResultFiles(self, resultDirs):
        mapping_result_files = []

        for result_dir in resultDirs:
            for root, dirs, files in os.walk(result_dir):
                for filename in files:
                    mapping_result_file = os.path.join(root, filename)
                    mapping_result_files.append(mapping_result_file)
        
        return mapping_result_files
                        

class BsfCall(BsfCallBase):
    """
    class to execute bsf-call process.
    """

    def __init__(self, refGenome, readFilePaths, cmdOpts):
        """
        constructor of BsfCall
        """

        self.refGenome = refGenome
        self.readFilePaths = readFilePaths
        self.reads = []
        self.opts = cmdOpts

        self.dataDir = None
        self.genomeDir = None
        self.mcContextDir = None

        self.readInFh1 = None
        self.readInFh2 = None
        self.numReads = {1: 0, 2: 0}

        self.mappingResultDirs = []
        self.mappingResultFiles = []

        self.setDataDir()
        self.setLogger()

        logging.info("bsf-call start.")
        self.logOption()

        self.numReadsPerFile = self.sizeForSplitRead(self.opts["split_read_size"])

        if self.opts["mapping_dir"]:
            self.opts["only_mcdetection"] = True
        else:
            self.opts["only_mcdetection"] = False


    def execute(self):
        """
        execute bsf-call process.
        """

        try:
            if self.opts["mapping_dir"]:
                # Only mc detection
                self.mappingResultDirs = self.opts["mapping_dir"].split(",")
                self.mappingResultFiles = self.getAllMappingResultFiles(self.mappingResultDirs)
                self.opts["mapping_result_files"] = self.mappingResultFiles
            else:
                self.makeIndexFile()
                self.processReads()
                self.mappingResultDirs = self.processMapping()

            self.processMcDetection(self.mappingResultDirs, self.opts["local_dir"])

            logging.info("bsf-call done.")
        except:
            logging.exception("Exception has occurred.")


    def processMapping(self):
        """
        run read mapping and filtering process.
        """

        logging.info("Mapping and filtering process start.")

        result_dirs = []
        for read_attr in self.reads:
            self.processMappingOneRead(read_attr)
            result_dirs.append(read_attr["results_dir"])

        logging.info("Mapping and filtering process done.")

        return result_dirs


    def processMcDetection(self, resultDirs, localDir = None):
        """
        run mC detection process.
        """

        logging.info("mC detection process start.")

        mc_detector = McDetector(self.refGenome, resultDirs, self.mcContextDir, self.opts)
        mc_detector.execute(self.opts["output"], self.opts["num_threads"])

        logging.info("mC detection process done.")


    def setDataDir(self):
        """
        create directries to store the files of the bsf-call process.
        """

        if self.opts["work_dir"]:
            self.dataDir = self.opts["work_dir"]
        else:
            self.dataDir = self.autoWorkDir()

        self.mcContextDir = "%s/mc_contexts" % self.dataDir

        if not os.path.exists(self.dataDir):
            os.makedirs(self.dataDir)

        if not os.path.exists(self.mcContextDir):
            os.makedirs(self.mcContextDir)


    def setLogger(self):
        """
        create logger to store the logs of the bsf-call process.
        """

        log_level = logging.INFO

        log_file = "%s/bsf-call.log" % self.dataDir
        file_logger = logging.FileHandler(filename=log_file)

        file_logger.setLevel(log_level)
        file_logger.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))

        logging.getLogger().addHandler(file_logger)
        logging.getLogger().setLevel(log_level)


    def logOption(self):
        """
        output bsf-call option values and arguments to the log.
        """

        if self.opts["mapping_dir"]:
            logging.info("Mapping result directory is specified. Only mC detection is executed.")
            logging.info("  Mapping result directory: %s" % self.opts["mapping_dir"])
            logging.info("  Reference genome: %s" % self.refGenome)
            logging.info("  Read BAM file: %s" % ("Yes" if self.opts["read_bam"] else "No"))
            logging.info("  Read SAM file: %s" % ("Yes" if self.opts["read_sam"] else "No"))
        else:
            logging.info("Reference genome: %s" % self.refGenome)
            logging.info("Read files: %s" % self.readFilePaths)
            logging.info("Working directory: %s" % self.dataDir)
            logging.info("Options:")
            logging.info("  Threshold of the alignment score at filtering: %d" % self.opts["aln_score_thres"])
            logging.info("  Paired-end direction: %s" % self.opts["pe_direction"])
            logging.info("  Options for LAST: %s" % self.opts["last_opts"])

        logging.info("  Threshold of read coverate: %d" % self.opts["coverage"])
        logging.info("  Threshold of mC ratio: %s" % str(self.opts["lower_bound"]))
        logging.info("  Threshold of the mismap probability at filtering: %s" % str(self.opts["aln_mismap_prob_thres"]))
        logging.info("  Working directory: %s" % self.dataDir)
        logging.info("  Local directory: %s" % self.opts["local_dir"])
        logging.info("  Output file: %s" % (self.opts["output"] if self.opts["output"] else "(stdout)"))
        logging.info("  Use cluster: %s" % ("Yes" if self.opts["use_cluster"] else "No"))
        logging.info("  Queue name: %s" % self.opts["queue_list"])
        logging.info("  Number of threads: %d" % self.opts["num_threads"])
        logging.info("  Split read size: %s" % self.opts["split_read_size"])


    def sizeForSplitRead(self, size):
        """
        get split read size.
        """

        unit = size[-1:].upper()
        if unit == "K":
            split_size = int(size[0:-1]) * 1000
        elif unit == "M":
            split_size = int(size[0:-1]) * 1000 * 1000
        else:
            split_size = int(size)

        return split_size


    def processReads(self):
        """
        process all read files.
        """

        for read_no, read_path in enumerate(self.readFilePaths):
            read = self.processOneRead(read_path, read_no + 1)
            self.reads.append(read)


    def processOneRead(self, readPath, readNo):
        """
        process one read. do the following:
        create directories to store the splited read and result files.
        if read is SRA file, convert to FASTQ file using fastq-dump command.
        split read file.
        """

        logging.info("Process read file start: %d: %s" % (readNo, readPath))

        data_dir = "%s/%d"% (self.dataDir, readNo)
        read = {"base_dir": data_dir, "path": readPath, "reads_dir": data_dir + "/reads", "results_dir": data_dir + "/results"}

        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
            os.makedirs(read["reads_dir"])
            os.makedirs(read["results_dir"])

        pe_no = self.pairedEndReadNumbers()[0]
        for read_path in readPath.split(","):
            dir_name, file_name, base_name, ext = self.splitFilePath(read_path)
            file_type = self.checkReadFileType(read_path)
            read[pe_no] = {"name": base_name, "fname": file_name, "type": file_type, "path": read_path}
            pe_no += 1

        if self.isPairedEnd(read):
            for pe_no in self.pairedEndReadNumbers():
                if read[pe_no]["type"] == "sra":
                    self.sra2Fastq(read[pe_no]["path"], data_dir, False)
                    fastq_fpath = self.fastqDumpedFilePath(data_dir, read[pe_no]["name"])
                    if os.path.exists(fastq_fpath):
                        read[pe_no]["name"] = self.readNameByReadFile(fastq_fpath)
                        read[pe_no]["path"] = fastq_fpath
                        read[pe_no]["type"] = "fastq"
                        read[pe_no]["fname"] = os.path.basename(read[pe_no]["path"])
        else:
            if read[1]["type"] == "sra":
                self.sra2Fastq(read[1]["path"], data_dir, True)
                read_name = read[1]["name"]
                for pe_no in self.pairedEndReadNumbers():
                    fastq_fpath = self.fastqDumpedFilePath(data_dir, read_name, pe_no)
                    if os.path.exists(fastq_fpath):
                        if not pe_no in read:
                            read[pe_no] = {}
                        read[pe_no]["name"] = self.readNameByReadFile(fastq_fpath)
                        read[pe_no]["path"] = fastq_fpath
                        read[pe_no]["type"] = "fastq"
                        read[pe_no]["fname"] = os.path.basename(read[pe_no]["path"])

        is_paired_end = self.isPairedEnd(read)
        logging.info("Paired-end: %s" % is_paired_end)
        logging.info("  Forward: %s" % read[1]["path"])
        if is_paired_end:
            logging.info("  Reverse: %s" % read[2]["path"])

        self.splitRead(read)

        # gzip FASTQ file which converted from SRA file in background.
        for dumped_fastq_file in glob.glob("%s/*.fastq" % read["base_dir"]):
            self.gzipFile(dumped_fastq_file, False)

        logging.info("Process read file done: %d: %s" % (readNo, readPath))

        return read


    def processMappingOneRead(self, readAttr):
        """
        run mapping and filtering process for one read.
        """

        logging.info("Target read file: %s" % readAttr["path"])

        self.runLast(readAttr)


    def runLast(self, readAttr):
        """
        run LAST programs to map reads and filtering.
        """

        is_paired_end = self.isPairedEnd(readAttr)
        filter_option = self.filterOpts(self.opts["aln_mismap_prob_thres"], self.opts["aln_score_thres"], is_paired_end)

        if is_paired_end:
            last_executor = LastExecutorPairedEnd(self.refGenome, self.dataDir, readAttr["reads_dir"], readAttr["results_dir"])
        else:
            last_executor = LastExecutorSingle(self.refGenome, self.dataDir, readAttr["reads_dir"], readAttr["results_dir"])

        last_executor.execute(readAttr, self.opts["num_threads"], self.lastalOpts(self.opts["last_opts"]), self.mergeOpts(), filter_option)


    def sra2Fastq(self, sraFile, outputDir, splitFiles = True):
        """
        convert SRA file to FASTQ file.
        """

        logging.info("Convert SRA file to FASTQ start: %s" % sraFile)

        cmd = "fastq-dump -O %s " % outputDir
        if splitFiles:
            cmd += "--split-files "
        cmd += sraFile
        logging.info(cmd)
        p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        out = p.stdout
        p.wait()
        out_data = out.read()
        if len(out_data) > 0:
            logging.info(out_data)

        logging.info("Convert SRA file to FASTQ done: %s" % sraFile)


    def splitRead(self, readAttr):
        """
        split read file.
        """

        logging.info("Split read file start: %s" % readAttr["path"])

        is_paired_end = self.isPairedEnd(readAttr)

        self.readInFh1 = open(readAttr[1]["path"], "r")
        if is_paired_end:
            self.readInFh2 = open(readAttr[2]["path"], "r")

        out_fpath1 = self.splitedReadFilePath(readAttr["reads_dir"], 1, self.numReadsPerFile, 1, readAttr[1]["type"])
        while self.outputOneSplittedReadFile(self.readInFh1, out_fpath1, readAttr, 1):
            out_fpath1 = self.splitedReadFilePath(readAttr["reads_dir"], self.numReads[1] + 1, self.numReads[1] + self.numReadsPerFile, 1, readAttr[1]["type"])

        if is_paired_end:
            out_fpath2 = self.secondReadFilePathByFirstReadFilePath(out_fpath1, readAttr[2]["type"])
            self.outputOneSplittedReadFile(self.readInFh2, out_fpath2, readAttr, 2)
            
        self.readInFh1.close()
        if is_paired_end:
            self.readInFh2.close()

        logging.info("Split read file done: %s" % readAttr["path"])
        logging.info("%s: the number of reads (1st read): %d" % (readAttr[1]["path"], self.numReads[1]))
        if is_paired_end:
            logging.info("%s: the number of reads (2nd read): %d" % (readAttr[2]["path"], self.numReads[2]))

        self.readInFh1 = None
        self.readInFh2 = None
        self.numReads[1] = 0
        self.numReads[2] = 0


    def outputOneSplittedReadFile(self, fh, outFilePath, readAttr, pairedEndReadNo):
        """
        output one splitted read file.
        """

        if readAttr[pairedEndReadNo]["type"] == "fastq":
            return self.outputOneSplittedFastqFile(fh, outFilePath, readAttr, pairedEndReadNo)
        elif readAttr[pairedEndReadNo]["type"] == "fasta":
            return self.outputOneSplittedFastaFile(fh, outFilePath, readAttr, pairedEndReadNo)


    def outputOneSplittedFastqFile(self, fh, outFilePath, readAttr, pairedEndReadNo):
        """
        output one splitted FASTQ file.
        """

        num_lines_per_read = 4
        line_no = 1
        out_fh = open(outFilePath, "w")
        for line in fh:
            if line_no % num_lines_per_read == 2:
                line = line.upper().replace("C", "t")
            out_fh.write(line)

            if line_no % num_lines_per_read == 0:
                self.numReads[pairedEndReadNo] += 1

            if line_no % (num_lines_per_read * self.numReadsPerFile) == 0:
                out_fh.close()
                self.afterProcessSplitRead(outFilePath, readAttr)
                if self.isPairedEnd(readAttr) and pairedEndReadNo == 1:
                    out_fpath2 = self.secondReadFilePathByFirstReadFilePath(outFilePath, readAttr[2]["type"])
                    self.outputOneSplittedReadFile(self.readInFh2, out_fpath2, readAttr, 2)

                return line_no

            line_no += 1
                
        self.afterProcessSplitRead(outFilePath, readAttr)

        return None


    def outputOneSplittedFastaFile(self, fh, outFilePath, readAttr, pairedEndReadNo):
        """
        output one splitted FASTA file.
        """

        num_reads = 0
        line_no = 1
        out_fh = open(outFilePath, "w")
        for line in fh:
            if line[0:1] == ">":
                out_fh.write(line)
                num_reads += 1
                self.numReads[pairedEndReadNo] += 1
            else:
                line = line.upper().replace("C", "t")
                out_fh.write(line)

                if num_reads != 0 and num_reads % self.numReadsPerFile == 0:
                    out_fh.close()
                    self.afterProcessSplitRead(outFilePath, readAttr)
                    if self.isPairedEnd(readAttr) and pairedEndReadNo == 1:
                        out_fpath2 = self.secondReadFilePathByFirstReadFilePath(outFilePath, readAttr[2]["type"])
                        self.outputOneSplittedReadFile(self.readInFh2, out_fpath2, readAttr, 2)
                    return line_no

            line_no += 1
                
        self.afterProcessSplitRead(outFilePath, readAttr)

        return None

                
    def afterProcessSplitRead(self, readFile, readAttr = None):
        """
        this function is called after output splitted one read file.
        """

        self.gzipFile(readFile, True, True)


    def makeIndexFile(self):
        """
        create index file of reference genome.
        """

        directions = []
        if not os.path.exists("%s.f.prj" % self.refGenome):
            directions.append("f")
        if not os.path.exists("%s.r.prj" % self.refGenome):
            directions.append("r")

        if len(directions) > 0:
            logging.info("Make index file start.")
            last_executor = LastExecutor(self.refGenome, self.dataDir)
            last_executor.lastdb(directions, self.opts["num_threads"] > 1)
            logging.info("Make index file done.")


    def checkReadFileType(self, readFilePath):
        """
        get read file type.
        """

        name, ext = os.path.splitext(readFilePath)
        if ext == ".gz":
            name, ext = os.path.splitext(name)

        if len(ext) > 1:
            ext = ext[1:]

        file_type = None

        if ext == "sra" or "fastq" or "fasta":
            file_type = ext
        elif ext == "fa":
            file_type = "fasta"
        else:
            f = open(readFilePath, "r")
            first_char = f.read(1)
            if first_char == "@":
                file_type = "fastq"
            elif first_char == ">":
                file_type = "fasta"
            else:
                file_type = "sra"
            f.close()

        return file_type


    def splitedReadFilePath(self, outputDir, start, end, readDirection, ext):
        """
        get splitted read file path.
        """

        return "%s/%010d-%010d_%d.%s" % (outputDir, start, end, readDirection, ext)


    def fastqDumpedFilePath(self, outputDir, readName, readDirection = None):
        """
        get output file path for fastq-dump command.
        """

        path = "%s/%s" % (outputDir, readName)
        if readDirection:
            path += "_%d" % readDirection 

        return path + ".fastq"


    def autoWorkDir(self):
        """
        get working directory path determined automatically.
        """

        now = datetime.now()

        s = ",".join(self.readFilePaths) + self.refGenome
        for key, value in self.opts.items():
            s += ("%s:%s" % (key, value))
        h = hashlib.md5(s).hexdigest()

        return "%s-%06d-%s" % (now.strftime("%Y%m%d-%H%M%S"), now.microsecond, h[0:16])
        

    def bisulfite_f_seed(self):
        """
        output bisulfite_f.seed file.
        """

        filename = 'bisulfite_f.seed'
        if not os.path.exists(filename):
            file = open(filename, "w")
            file.write("""
# This subset seed pattern is suitable for aligning forward strands of
# bisulfite-converted DNA to unconverted DNA.
CT A G
CT A G
CT A G
CT A G
CT A G
CT A G
CTAG
CT A G
CTAG
CT A G
CT A G
CTAG
CTAG
""")
            file.close()
        return


    def bisulfite_r_seed(self):
        """
        output bisulfite_r.seed file.
        """

        filename = 'bisulfite_r.seed'
        if not os.path.exists(filename):
            file = open(filename, "w")
            file.write("""
# This subset seed pattern is suitable for aligning reverse strands of
# bisulfite-converted DNA to unconverted DNA.
AG C T
AG C T
AG C T
AG C T
AG C T
AG C T
AGCT
AG C T
AGCT
AG C T
AG C T
AGCT
AGCT
""")
            file.close()
        return


    def bisulfite_f_mat(self):
        """
        output bisulfite_f.mat file.
        """

        filename = 'bisulfite_f.mat'
        if not os.path.exists(filename):
            file = open(filename, "w")
            file.write("""
# This score matrix is suitable for aligning forward strands of
# bisulfite-converted DNA queries to unconverted reference sequences.
    A   C   G   T
A   6 -18 -18 -18
C -18   6 -18   3
G -18 -18   6 -18
T -18 -18 -18   3
""")
            file.close()
        return


    def bisulfite_r_mat(self):
        """
        output bisulfite_r.mat file.
        """

        filename = 'bisulfite_r.mat'
        if not os.path.exists(filename):
            file = open(filename, "w")
            file.write("""
# This score matrix is suitable for aligning forward strands of
# bisulfite-converted DNA queries to unconverted reference sequences.
    A   C   G   T
A   3 -18 -18 -18
C -18   6 -18 -18
G   3 -18   6 -18
T -18 -18 -18   6
""")
            file.close()
        return


class BsfCallCluster(BsfCall):
    """
    class to execute bsf-call process on pc cluster.
    """

    def __init__(self, refGenome, readFilePaths, cmdOpts):
        """
        constructor of BsfCallCluster
        """

        BsfCall.__init__(self, refGenome, readFilePaths, cmdOpts)

        self.lastExecutor = LastExecutorCluster(self.refGenome, self.opts)
        self.mappingJobIds = []


    def processMapping(self):
        """
        if bsf-call is executed on cluster, mapping and filtering job is submitted
        when read file is splitted.
        therefore, in this function, only wait for all jobs to finish.
        """

        logging.info("Waiting for all jobs to finish.")

        self.waitForSubmitJobs(self.mappingJobIds)

        logging.info("Mapping and filtering process done.")

        return self.lastExecutor.resultDirs


    def processMcDetection(self, resultDirs, localDir = None):
        """
        for each chromosome number, submit mC detection process job to the cluster.
        after all jobs have been finished, output mC detection result.
        """

        logging.info("mC detection process start.")

        chrs = self.chrsByRefGenome(self.refGenome)
        job_ids = []
        for chr_no in chrs:
            job_id = self.submitMcDetectionJob(resultDirs, chr_no)
            job_ids.append(job_id)
            sleep(1)
        logging.info("Submitted jobs: %s" % ",".join(job_ids))

        self.waitForSubmitJobs(job_ids)

        mc_detector = McDetector(self.refGenome, resultDirs, self.mcContextDir, self.opts)        

        mc_detector.chrs = chrs
        mc_detector.output(self.opts["output"])

        logging.info("mC detection process done.")


    def submitMcDetectionJob(self, resultDirs, chrNo):
        """
        submit mC detection process job to the cluster.
        """

        argv = self.qsubRemoteCommandArgv(resultDirs, chrNo)
        cmd = self.qsubCommand(chrNo, " ".join(map((lambda s: '"' + s + '"'), argv)))

        qsub = os.popen(cmd)
        out = qsub.read()
        qsub.close()

        job_id = self.jobIdByQsubOutput(out)
        self.logJobSubmit("mC detection job: %s: Mapping result directories: %s" % (chrNo, resultDirs), job_id)

        return job_id


    def qsubRemoteCommandArgv(self, resultDirs, chrNo):
        """
        get arguments for mC detection program (mc-detector.py).
        """

        argv = []

        argv.append(self.refGenome)
        argv.append(",".join(resultDirs))
        argv.append(self.mcContextDir)
        argv.append(chrNo)
        argv.append(str(self.opts["only_mcdetection"]))
        argv.append(str(self.opts["lower_bound"]))
        argv.append(str(self.opts["coverage"]))
        argv.append(str(self.opts["aln_mismap_prob_thres"]))

        if self.opts["read_bam"]:
            argv.append("BAM")
        elif self.opts["read_sam"]:
            argv.append("SAM")
        else:
            argv.append("MAF")
            
        if self.opts["local_dir"]:
            argv.append(self.opts["local_dir"])

        return argv


    def qsubCommand(self, chrNo, cmdArgs):
        """
        get qsub command to submit mC detection job.
        """

        remote_cmd = os.path.join(self.scriptDir(), "mc-detector.py")
        out_file = os.path.join(self.mcContextDir, "mc-detector-%s.out" % chrNo)
        err_file = os.path.join(self.mcContextDir, "mc-detector-%s.err" % chrNo)

        if self.opts["queue_list"]:
            cmd = "qsub -o %s -e %s -q %s -cwd -b y %s %s" % (out_file, err_file, self.opts["queue_list"], remote_cmd, cmdArgs)
        else:
            cmd = "qsub -o %s -e %s -cwd -b y %s %s" % (out_file, err_file, remote_cmd, cmdArgs)

        return cmd


    def afterProcessSplitRead(self, readFile, readAttr = None):
        """
        this function is called after output splitted one read file.
        if bsf-call is executed on pc cluster, mapping and filtering job is submitted
        after output one splitted read file.
        """

        if self.isPairedEnd(readAttr):
            fpath, ext = os.path.splitext(readFile)
            if fpath[-1:] == "2":
                read_file = "%s1.%s" % (fpath[0:-1], readAttr[1]["type"])
                job_id = self.lastExecutor.executeOneRead(readFile, readAttr)
                self.mappingJobIds.append(job_id)
        else:
            job_id = self.lastExecutor.executeOneRead(readFile, readAttr)
            self.mappingJobIds.append(job_id)


class LastExecutor(BsfCallBase):
    """
    class to run LAST programs to map read and filtering.
    """

    def __init__(self, refGenome, baseDir = ".", readsDir = None, resultsDir = None):
        """
        constructor of LastExecutor
        """

        self.refGenome = refGenome
        self.baseDir = baseDir
        self.readsDir = readsDir
        self.resultsDir = resultsDir
        self.queue = Queue.Queue()
        self.lock = threading.Lock()


    def execute(self, readAttr, numThreads, lastalOpts = "", mergeOpts = "", filterOpts = ""):
        """
        enqueue all splited read files to the queue.
        create and start threads to execute read mapping and filtering process.
        wait for all threads to finish.
        """

        self.enqueue(readAttr)

        threads = []
        for i in range(numThreads):
            t = threading.Thread(target=self.worker, args=(readAttr, lastalOpts, mergeOpts, filterOpts))
            t.daemon = True
            threads.append(t)
            t.start()

        for thread in threads:
            thread.join()


    def worker(self, readAttr, lastalOpts, mergeOpts, filterOpts):
        """
        thread worker to execute LAST.
        dequeue read file path from the queue and execute read mapping filtering
        process.
        """

        while True:
            try:
                fpath = self.queue.get_nowait()
                self.runLast(fpath, readAttr, lastalOpts, mergeOpts, filterOpts)
            except Queue.Empty:
                break

        
    def runLast(self, readFile, readAttr, lastalOpts, mergeOpts, filterOpts, rmInFiles = True):
        """
        execute LAST programs to map read and filtering.
        """

        cmd = self.batchCmd(readFile, readAttr, lastalOpts, mergeOpts, filterOpts, rmInFiles)
        logging.debug(cmd)
        p = subprocess.Popen(cmd, shell = True, stderr = subprocess.PIPE)
        error = p.stderr
        p.wait()

        error_msg = error.read()
        if len(error_msg) > 0:
            logging.fatal(error_msg)


    def enqueue(self, readAttr):
        """
        enqueue all splitted read files to the queue.
        """

        for read_file in glob.glob("%s/*_1.%s.gz" % (self.readsDir, readAttr[1]["type"])):
            self.queue.put(read_file)


    def lastdb(self, directions, parallel = False):
        """
        execute lastdb command to create index file of reference genome.
        """

        cmds = []
        for direction in directions:
            cmds.append(self.lastdbCmd(direction))

        if parallel:
            processes = []
            for cmd in cmds:
                logging.info(cmd)
                p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
                out = p.stdout
                processes.append({"p": p, "out": out})
            for process in processes:
                process["p"].wait()
                out_data = process["out"].read()
                if len(out_data) > 0:
                    logging.info(out_data)
        else:
            for cmd in cmds:
                logging.info(cmd)
                p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
                out = p.stdout
                p.wait()
                out_data = out.read()
                if len(out_data) > 0:
                    logging.info(out_data)


    def lastdbCmd(self, direction):
        """
        get lastdb command to create index file of reference genome.
        """

        return "lastdb -w2 -u bisulfite_%s.seed %s.%s %s" % (direction, self.refGenome, direction, self.refGenome)


    def lastalCmd(self, readFile, direction, opts = ""):
        """
        get lastal command to map read.
        """

        s_opt = self.lastalSopt(direction)
        read_name = self.readNameByReadFile(readFile)

        return "zcat %s | lastal -p bisulfite_%s.mat %s %s %s.%s - > %s" % (readFile, direction, s_opt, opts, self.refGenome, direction, self.mappingResultFilePath(read_name, direction))


    def mergeCmd(self, forwardFile, reverseFile, outputFile, opts = "", rmInFiles = True):
        """
        get command to merge lastal output.
        """

        cmd = "last-merge-batches %s %s %s > %s" % (opts, forwardFile, reverseFile, outputFile)
        if rmInFiles:
            cmd += "; rm %s %s" % (forwardFile, reverseFile)

        return cmd


    def mappingAndMergeCmd(self, readFile, lastalOpts = "", mergeOpts = "", rmInFiles = True):
        """
        get read mapping and filtering command.
        """

        read_name = self.readNameByReadFile(readFile)
        n, ext = os.path.splitext(readFile)
        if ext == ".gz":
            n, ext = os.path.splitext(n)
        lastal_qopt = self.lastalQopt(ext[1:])
        lastal_opt = "%s %s" % (lastalOpts, lastal_qopt)
        mapping_file_f = self.mappingResultFilePath(read_name, "f")
        mapping_file_r = self.mappingResultFilePath(read_name, "r")
        merged_file = self.mergeResultFilePath(read_name)

        return "%s; %s; %s" % (self.lastalCmd(readFile, "f", lastal_opt), self.lastalCmd(readFile, "r", lastal_opt), self.mergeCmd(mapping_file_f, mapping_file_r, merged_file, mergeOpts, rmInFiles))


    def mappingResultFilePath(self, readName, direction):
        """
        get read mapping result file path.
        """

        return "%s/%s_%s" % (self.resultsDir, readName, direction)


    def mergeResultFilePath(self, readName):
        """
        get merge result file path.
        """

        return "%s/%s" % (self.resultsDir, readName)


    def filterResultFilePath(self, readName):
        """
        get filtering result file path.
        """

        return "%s/%s.maf" % (self.resultsDir, readName)


    def lastalSopt(self, direction):
        """
        get -s option for lastal.
        """

        opt = ""
        if direction == "f":
            opt = "-s1"
        elif direction == "r":
            opt = "-s0"

        return opt


    def lastalQopt(self, fileType):
        """
        get -Q option for lastal.
        """

        opt = ""
        if fileType == "fasta":
            opt = "-Q0"
        elif fileType == "fastq":
            opt = "-Q1"

        return opt


class LastExecutorCluster(LastExecutor):
    """
    class to run LAST programs on pc cluster.
    """

    def __init__(self, refGenome, bsfCallOpts):
        """
        constructor of LastExecutorCluster
        """

        self.refGenome = refGenome
        self.opts = bsfCallOpts

        self.resultDirs = []


    def executeOneRead(self, readFile, readAttr):
        """
        execute read mapping and filtering process with specified read.
        on pc cluster, submit mapping and filtering job and return job id.
        """

        if readAttr["results_dir"] not in self.resultDirs:
            self.resultDirs.append(readAttr["results_dir"])

        lastal_opts = self.lastalOpts(self.opts["last_opts"])
        merge_opts = self.mergeOpts()
        filter_opts = self.filterOpts(self.opts["aln_mismap_prob_thres"], self.opts["aln_score_thres"], self.isPairedEnd(readAttr))
        job_id = self.submitJob(readFile, readAttr, lastal_opts, merge_opts, filter_opts, self.opts["queue_list"])

        return job_id


    def submitJob(self, readFile, readAttr, lastalOpts = "", mergeOpts = "", filterOpts = "", queueName = None):
        """
        submit read mapping and filtering process job.
        """

        job_id = None

        read_name = self.readNameByReadFile(readFile)[0:-2]
        out_file = self.qsubStdoutFilePath(readAttr["base_dir"], read_name)
        err_file = self.qsubStderrFilePath(readAttr["base_dir"], read_name)
        remote_cmd = self.remoteCommand(readAttr)
        remote_cmd_args = " ".join(map((lambda s: '"' + s + '"'), self.remoteCommandArgv(read_name, readAttr, lastalOpts, filterOpts)))

        if queueName:
            cmd = "qsub -o %s -e %s -q %s -cwd %s %s" % (out_file, err_file, queueName, remote_cmd, remote_cmd_args)
        else:
            cmd = "qsub -o %s -e %s -cwd %s %s" % (out_file, err_file, remote_cmd, remote_cmd_args)

        qsub = os.popen(cmd)
        out = qsub.read()
        qsub.close()

        job_id = self.jobIdByQsubOutput(out)

        dir_path, file_name, base_name, ext = self.splitFilePath(readFile)
        self.logJobSubmit("Mapping and filtering: read: %s/%s" % (dir_path, base_name[0:-2]), job_id)

        return job_id


    def remoteCommand(self, readAttr):
        """
        get read mapping and filtering command path to submit by qsub command.
        """

        if self.isPairedEnd(readAttr):
            return os.path.join(self.scriptDir(), "mapping-p.sh")           
        else:
            return os.path.join(self.scriptDir(), "mapping-s.sh")


    def remoteCommandArgv(self, readName, readAttr, lastalOpts, filterOpts):
        """
        get read mapping and filtering command arguments.
        """

        argv = []

        argv.append(readAttr["reads_dir"])
        argv.append(readAttr["results_dir"])
        argv.append(self.refGenome)
        argv.append(filterOpts)

        argv.append(readName)
        argv.append(readAttr[1]["type"])
        argv.append("%s %s" % (lastalOpts, self.lastalQopt(readAttr[1]["type"])))

        if self.isPairedEnd(readAttr):
            argv.append(readName)
            argv.append(readAttr[2]["type"])
            argv.append("%s %s" % (lastalOpts, self.lastalQopt(readAttr[2]["type"])))

        return argv


    def qsubStdoutFilePath(self, dirPath, readName):
        """
        get qsub command standard output file path.
        """

        return "%s/mapping_%s.out" % (dirPath, readName)


    def qsubStderrFilePath(self, dirPath, readName):
        """
        get qsub command standard error file path.
        """

        return "%s/mapping_%s.err" % (dirPath, readName)


class LastExecutorSingle(LastExecutor):
    """
    class to run LAST programs to map single read and filtering.
    """

    def __init__(self, refGenome, baseDir, readsDir, resultsDir):
        """
        constructor of LastExecutorSingle
        """

        LastExecutor.__init__(self, refGenome, baseDir, readsDir, resultsDir)


    def batchCmd(self, readFile, readAttr, lastalOpts = "", mergeOpts = "", filterOpts = "", rmInFiles = True):
        """
        get batch command to map read and filtering.
        """

        read_name = self.readNameByReadFile(readFile)
        out_file = self.filterResultFilePath(read_name[0:-2])

        cmds = []
        cmds.append(self.mappingAndMergeCmd(readFile, lastalOpts, mergeOpts, rmInFiles))
        cmds.append(self.filterCmd(self.mergeResultFilePath(read_name), out_file, filterOpts, rmInFiles))
        cmds.append("gzip %s" % out_file)

        return "; ".join(cmds)


    def filterCmd(self, inputFile, outputFile, opts = "", rmInFile = True):
        """
        get filter command.
        """

        cmd = "last-map-probs %s %s > %s" % (opts, inputFile, outputFile)
        if rmInFile:
            cmd += "; rm %s" % inputFile

        return cmd


class LastExecutorPairedEnd(LastExecutor):
    """
    class to run LAST programs to map paired-end read and filtering.
    """

    def __init__(self, refGenome, baseDir, readsDir, resultsDir):
        """
        constructor of LastExecutorPairedEnd
        """

        LastExecutor.__init__(self, refGenome, baseDir, readsDir, resultsDir)


    def batchCmd(self, readFile1, readAttr, lastalOpts = "", mergeOpts = "", filterOpts = "", rmInFiles = True):
        """
        get batch command to map read and filtering.
        """

        read_file2 = self.secondReadFilePathByFirstReadFilePath(readFile1, readAttr[2]["type"])
        read_files = (readFile1, read_file2)

        filter_cmd_in_file = "%s %s" % (self.mergeResultFilePath(read_files[0]), self.mergeResultFilePath(read_files[1]))
        read_name1 = self.readNameByReadFile(read_files[0])
        read_name2 = self.readNameByReadFile(read_files[1])
        merge_result_file = "%s %s" % (self.mergeResultFilePath(read_name1), self.mergeResultFilePath(read_name2))
        out_file = self.filterResultFilePath(read_name1[0:-2])

        cmds = []
        cmds.append(self.mappingAndMergeCmd(read_files[0], lastalOpts, mergeOpts, rmInFiles))
        cmds.append(self.mappingAndMergeCmd(read_files[1], lastalOpts, mergeOpts, rmInFiles))
        cmds.append(self.filterCmd(merge_result_file, out_file, filterOpts, rmInFiles))
        cmds.append("gzip %s" % out_file)

        return "; ".join(cmds)


    def filterCmd(self, inputFile, outputFile, opts = "", rmInFile = True):
        """
        get filter command.
        """

        cmd = "last-pair-probs %s %s > %s" % (opts, inputFile, outputFile)
        if rmInFile:
            cmd += "; rm %s" % inputFile

        return cmd


class McDetector(BsfCallBase):
    """
    class to execute mC detection process.
    """

    def __init__(self, refGenome, resultDirs, mcContextDir, options):
        """
        constructor of McDetector
        """

        self.refGenome = refGenome
        self.mappingResultDirs = resultDirs
        self.mcContextDir = mcContextDir
        self.lowerBound = options["lower_bound"]
        self.coverageThreshold = options["coverage"]
        self.onlyMcDetection = options["only_mcdetection"]

        self.opts = options

        self.mappingResultFiles = []

        if self.onlyMcDetection:
            self.mismapThreshold = options["aln_mismap_prob_thres"]
            self.readBam = options["read_bam"]
            self.readSam = options["read_sam"]
            if "mapping_result_files" in options:
                self.mappingResultFiles = options["mapping_result_files"]
            else:
                self.mappingResultFiles = self.getAllMappingResultFiles(resultDirs)
        else:
            self.mismapThreshold = None
            self.readBam = False
            self.readSam = False

        if options["local_dir"]:
            self.localDir = options["local_dir"]
        else:
            self.localDir = mcContextDir

        self.chrs = []
        self.refGenomeBuf = []
        self.targetChr = ""
        self.targetSeqD = ""
        self.targetSeqC = ""
        self.targetSeqLen = 0

        self.mcDetectData = {}

        logging.debug("McDetector.__init__: self.mappingResultDirs: %s" % self.mappingResultDirs)


    def execute(self, outputFile, numWorkers = 1):
        """
        execute mC detection process and output result.
        """

        for line in open(self.refGenome, "r"):
            line = line.strip()
            if line[0] == ">":
                if self.targetChr != "":
                    self.processOneChr()
                    self.chrs.append(self.targetChr)
                self.targetChr = self.chrnoFromFastaDescription(line)
            else:
                self.refGenomeBuf.append(line.upper())
        # last chromosome
        self.processOneChr()
        self.chrs.append(self.targetChr)

        self.output(outputFile)


    def processOneChr(self, chrNo = None):
        """
        execute mC detection process with specified chromosome number.
        """

        if chrNo:
            self.targetChr = chrNo

        local_dir = "%s/%s" % (self.localDir, self.targetChr)
        if not os.path.exists(local_dir):
            os.makedirs(local_dir)

        if self.refGenomeBuf:
            self.targetSeqD = "".join(self.refGenomeBuf)
            del self.refGenomeBuf[:]
        else:
            self.targetSeqD = self.oneChrGenomeSeq(self.refGenome, chrNo)

        self.targetSeqD = self.targetSeqD.upper()
        self.targetSeqLen = len(self.targetSeqD)
        self.targetSeqC = self.complementSeq(self.targetSeqD).upper()

        logging.info("mC detection process start: %s (%d)" % (self.targetChr, self.targetSeqLen))

        for mapping_result_file in self.mappingResultFiles:
            logging.debug("Mapping result file: %s" % mapping_result_file)
            if self.onlyMcDetection:
                if not (self.readBam or self.readSam):
                    if self.isGzippedFile(mapping_result_file) or self.isMafFile(mapping_result_file):
                        self.processMafFile(mapping_result_file)
                else:
                    # BAM or SAM
                    if self.readBam and self.isBamFile(mapping_result_file):
                        self.processSamFile(mapping_result_file)
                    elif self.readSam and self.isSamFile(mapping_result_file):
                        self.processSamFile(mapping_result_file)
            else:
                self.processMafFile(mapping_result_file)

        logging.info("mC detection process done: %s (%d)" % (self.targetChr, self.targetSeqLen))

        logging.info("Output result file start.")
        self.outputOneChr(self.targetChr)
        logging.info("Output result file done.")

        if self.mcContextDir != self.localDir:
            copy(self.mcCotextLocalFilePath(self.targetChr), self.mcContextDir)

        self.targetSeqD = ""
        self.targetSeqC = ""
        self.targetSeqLen = 0


    def output(self, outputFile):
        """
        output mC detection result.
        """

        popen_args = ['cat']
        for chr in sorted(self.chrs, self.chrSort):
            result_file = self.mcContextGlobalFilePath(chr)
            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                popen_args.append(result_file)

        if len(popen_args) == 1:
            if outputFile:
                f = open(outputFile, "w")
                f.close()
        else:
            if outputFile:
                output_fh = open(outputFile, "w")
            else:
                output_fh = None
            subprocess.check_call(popen_args, stdout = output_fh)
            if output_fh:
                output_fh.close()


    def outputOneChr(self, chrNo):
        """
        output mC detection result with specified chromosome number.
        """

        logging.debug("McDetector.outputOneChr: chr = %s" % chrNo)
        ctx_dir = "%s/%s" % (self.localDir, chrNo)
        if os.path.isdir(ctx_dir):
            for fname in sorted(os.listdir(ctx_dir)):
                self.outputByOneMcContextFile(chrNo, "%s/%s" % (ctx_dir, fname))


    def outputByOneMcContextFile(self, chrNo, fpath):
        """
        output mC detection result with specified mC context file.
        """

        mc_contexts = {}
        f = open(fpath, "r")
        for line in f:
            try:
                ctx_pos, strand, mc_ctx, base = line.strip().split("\t")
                ctx_pos = int(ctx_pos)
                if not ctx_pos in mc_contexts:
                    mc_contexts[ctx_pos] = {}
                if not strand in mc_contexts[ctx_pos]:
                    mc_contexts[ctx_pos][strand] = {}
                if not mc_ctx in mc_contexts[ctx_pos][strand]:
                    mc_contexts[ctx_pos][strand][mc_ctx] = []
                mc_contexts[ctx_pos][strand][mc_ctx].append(base)
            except ValueError, e:
                logging.warning("ValueError: %s: %s -> %s" % (fpath, line.strip(), e.args[0]))
        f.close()

        out_f = open(self.mcCotextLocalFilePath(chrNo), "a")
        for pos in sorted(mc_contexts.keys()):
            for strand in mc_contexts[pos].keys():
                for mc_ctx in mc_contexts[pos][strand].keys():
                    coverage, mc_ratio = self.calcCoverage(mc_contexts[pos][strand][mc_ctx])
                    if coverage >= self.coverageThreshold and mc_ratio >= self.lowerBound:
                        out_f.write("%s\t%d\t%s\t%s\t%.02f\t%d\n" % (chrNo, pos, strand, mc_ctx, mc_ratio, coverage))
        out_f.close()

        self.gzipFile(fpath, False)


    def processMafFile(self, resultFile):
        """
        read mapping result file, and output mC context data file.
        """

        logging.info("Process MAF file start: %s" % resultFile)

        chr_line = True
        chr = ""
        f = None
        has_mismap = True
        mismap_value = None
        skip = False
        try:
            if self.isGzippedFile(resultFile):
                f = gzip.open(resultFile, "rb")
            else:
                f = open(resultFile, "r")

            counter = 1
            for line in f:
                tag = line[0]
                if tag != "s" and tag != "a":
                    continue

                if tag == "s" and skip:
                    continue

                line = line.strip()

                if tag == "a":
                    if self.mismapThreshold is not None:
                        skip = False
                        mismap_value = self.mafMismapValue(line)
                        if mismap_value is None:
                            logging.debug("mismap probability: None")
                            has_mismap = False
                        else:
                            logging.debug("mismap probability: %f" % mismap_value)
                            if mismap_value - self.mismapThreshold > 0:
                                skip = True
                    continue

                tag, name, start, aln_size, strand, seq_size, alignment = line.split()
                if chr_line:
                    if name != self.targetChr:
                        continue
                    chr = name
                    gl_alignment = self.clearGap(alignment)
                    ref_subseqlen = len(gl_alignment)
                    ref_start = int(start)
                    ref_seqlen = int(seq_size)
                    ref_seq = alignment
                else:
                    if chr == "":
                        continue
                    if strand == "-":
                        ref_start = self.complementStartPosition(ref_seqlen, ref_start, ref_subseqlen)
                        ref_seq = self.complementSeq(ref_seq)
                        alignment = self.complementSeq(alignment)
                    read_seq = alignment.replace("t", "C")
                    self.extractMcContextsByOneRead(chr, strand, ref_seq.upper(), ref_start, ref_seqlen, read_seq)
                    counter += 1
                    if counter > 500000:
                        self.outputMcContextData()
                        counter = 1

                chr_line = not chr_line
            f.close()

            self.outputMcContextData()

            if self.mismapThreshold and not has_mismap:
                logging.warning("MAF file has no mis-mapping probability.")
            logging.info("Process MAF file done: %s" % resultFile)
        except Exception, e:
            logging.fatal(resultFile)
            logging.fatal("  %s" % str(type(e)))
            logging.fatal("  %s" % str(e))
            if f:
                f.close()


    def processSamFile(self, samFile):
        """
        read mapping BAM/SAM file, and output mC context data file.
        """

        logging.info("Process BAM/SAM file start: %s" % samFile)

        samfile = None
        try:
            if self.readBam:
                samfile = pysam.Samfile(samFile, "rb")
            else:
                samfile = pysam.Samfile(samFile, "r")

            counter = 1
            for aln in samfile.fetch(until_eof = True):
                samaln = SamAlnParser(samfile, aln)
                samaln.setRefGenome(self.targetChr, self.targetSeqD, self.targetSeqLen)
                samaln.parse()

                chr = samaln.referenceName()

                if (chr is None) or (chr != self.targetChr):
                    continue

                if aln.mapq:
                    mismap = self.bamMapq2Mismap(aln.mapq)
                    if mismap - self.mismapThreshold > 0:
                        continue

                read_seq = samaln.alnReadSeq.replace("t", "C")
                self.extractMcContextsByOneRead(chr, samaln.strand, samaln.alnRefSeq.upper(), samaln.alnRefStart, self.targetSeqLen, read_seq)

                counter += 1
                if counter > 500000:
                    self.outputMcContextData()
                    counter = 1

            samfile.close()
            self.outputMcContextData()
        except Exception, e:
            logging.fatal(samFile)
            logging.fatal("  %s" % str(type(e)))
            logging.fatal("  %s" % str(e))
            if samfile:
                samfile.close()


    def extractMcContextsByOneRead(self, chr, strand, refSeq, refStart, refSeqLen, readSeq):
        """
        extract mC context by one read.
        """

        logging.debug("extractMcContextsByOneRead(%s, %s, %s, %d, %d, %s)" % (chr, strand, refSeq, refStart, refSeqLen, readSeq))

        nogap_refseq = self.clearGap(refSeq)
        bases = list(refSeq)
        last_pos = len(nogap_refseq) - 1
        pos = -1
        while True:
            try:
                pos = bases.index("C", pos + 1)
                num_gaps = refSeq.count("-", 0, pos)
                real_pos = pos - num_gaps
                ctx_type = self.mcContextType(nogap_refseq, real_pos)
                ctx_pos = refStart + real_pos
                if ctx_type == None:
                    if strand == "+":
                        ctx_type = self.mcContextType(self.targetSeqD, ctx_pos)
                    elif strand == "-":
                        ctx_type = self.mcContextType(self.targetSeqC, ctx_pos)
                if strand == "-":
                    ctx_pos = refSeqLen - ctx_pos - 1
                line = "%d\t%s\t%s\t%s\n" % (ctx_pos, strand, ctx_type, readSeq[pos])
                base_name = self.mcContextFileBaseName(ctx_pos)
                if base_name not in self.mcDetectData:
                    self.mcDetectData[base_name] = []
                self.mcDetectData[base_name].append(line)
            except IndexError:
                logging.debug("extractMcContextsByOneRead#IndexError: %s %d %s %s %s %d" % (chr, ctx_pos, strand, ctx_type, readSeq, pos))
            except ValueError:
                break


    def outputMcContextData(self):
        """
        output mC context data to file.
        """

        logging.info("Output mC detection data start.")

        for base_name in self.mcDetectData.keys():
            fpath = self.mcContextFilePathByName(base_name)
            logging.debug("%s: %d" % (fpath, len(self.mcDetectData[base_name])))
            fo = open(fpath, "a")
            fo.write("".join(self.mcDetectData[base_name]))
            fo.close()
        self.mcDetectData.clear()

        logging.info("Output mC detection data done.")


    def mcContextHash(self):
        h = {}
        for strand in self.strands():
            h[strand] = {}
            for mc_ctx in self.mcContextTypes():
                h[strand][mc_ctx] = {}

        return h


    def isC(self, seq):
        return seq[0:1].upper() == "C"


    def isT(self, seq):
        return seq[0:1].upper() == "T"


    def calcCoverage(self, seqAry):
        """
        count the number of C and T, calculate mC ratio.
        """

        num_seq = len(seqAry)
        c_ary = filter(self.isC, seqAry)
        t_ary = filter(self.isT, seqAry)
        num_c = len(c_ary)
        num_t = len(t_ary)

        if num_c + num_t == 0:
            return (num_seq, 0)
        else:
            return (num_seq, float(num_c) / (float(num_c) + float(num_t)))


    def mafMismapValue(self, aLine):
        """
        get mismap value by maf "a" line.
        """

        mismap_fields = filter((lambda s: s.startswith('mismap=')), aLine.split())
        if len(mismap_fields) > 0:
            return float(mismap_fields[0][7:])
        else:
            return None


    def mcContextFilePath(self, pos):
        """
        get mC context data file path with specified position.
        """

        return self.mcContextFilePathByName(self.mcContextFileBaseName(pos))


    def mcContextFilePathByName(self, name):
        """
        get mC context data file path with specified name.
        """

        return "%s/%s/%s.tsv" % (self.localDir, self.targetChr, name)


    def mcContextFileBaseName(self, pos):
        """
        get mC context data file name with specified chromosome number.
        """

        return "%010d" % (int(pos) / 100000)


    def mcCotextLocalFilePath(self, chrNo):
        """
        get local mC context data file path with specified chromosome number.
        """

        return "%s/%s.tsv" % (self.localDir, chrNo)


    def mcContextGlobalFilePath(self, chrNo):
        """
        get global mC context data file path with specified chromosome number.
        """

        return "%s/%s.tsv" % (self.mcContextDir, chrNo)


class SamAlnParser(BsfCallBase):

    def __init__(self, samfile, aln):
        self.samfile = samfile
        self.aln = aln

        self.refName = None
        self.refSeq = None
        self.refSeqLen = None

        self.strand = None

        self.alnRefStart = None
        self.alnRefSeq = None
        self.alnRefSeqLen = None

        self.alnReadSeq = None


    def setRefGenome(self, name, sequence, length = None):
        if length is None:
            length = len(sequence)

        self.refName = name
        self.refSeq = sequence
        self.refSeqLen = length


    def referenceName(self):
        if self.aln.tid >= 0:
            return self.samfile.getrname(self.aln.tid)
        else:
            return None


    def getStrand(self):
        if self.aln.is_reverse:
            return "-"
        else:
            return "+"


    def parseCigar(self, cigar):
        return [{"op": v[-1:], "num": int(v[:-1])} for v in re.findall('\d+[A-Z=]', cigar)]


    def alignmentSequences(self, refSeqPos, readSeq, cigars):
        refs = []
        rest_refseq = 0

        reads = []
        read_pos = 0

        for cigar in cigars:
            if cigar["op"] == "M":
                refs.append(self.refSeq[refSeqPos:refSeqPos+cigar["num"]])
                refSeqPos += cigar["num"]
                reads.append(readSeq[read_pos:read_pos+cigar["num"]])
                read_pos += cigar["num"]
            elif cigar["op"] == "I":
                refs.append("-" * cigar["num"])
                reads.append(readSeq[read_pos:read_pos+cigar["num"]])
                read_pos += cigar["num"]
            elif cigar["op"] == "P":
                refs.append("-" * cigar["num"])
                reads.append("-" * cigar["num"])
            elif cigar["op"] == "D" or cigar["op"] == "N":
                rest_refseq += cigar["num"]
                reads.append("-" * cigar["num"])
            elif cigar["op"] == "S":
                read_pos += cigar["num"]

        if rest_refseq > 0:
            refs.append(self.refSeq[refSeqPos:refSeqPos+rest_refseq])

        return {"reference": "".join(refs), "read": "".join(reads)}


    def parse(self):
        self.strand = self.getStrand()

        self.alnRefStart = self.aln.pos
        read_seq = self.aln.seq
        cigars = self.parseCigar(self.aln.cigarstring)
        alignment = self.alignmentSequences(self.alnRefStart, read_seq, cigars)

        if self.strand == "+":
            self.alnRefSeq = alignment["reference"]
            self.alnReadSeq = alignment["read"]
        elif self.strand == "-":
            nogap_refseq = self.clearGap(alignment["reference"])
            self.alnRefStart =  self.complementStartPosition(self.refSeqLen, self.alnRefStart, len(nogap_refseq))
            self.alnRefSeq = self.complementSeq(alignment["reference"])
            self.alnReadSeq = self.complementSeq(alignment["read"])

