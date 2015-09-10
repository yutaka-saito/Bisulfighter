#!/usr/bin/env python
"""
Bisulfighter::bsf-call

Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
by National Institute of Advanced Industrial Science and Technology (AIST)
is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
http://creativecommons.org/licenses/by-nc-sa/3.0/
"""

__version__= "1.3"

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
import bz2
import zipfile
import logging
from time import sleep
from shutil import copy
import re
# import pysam

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
        dir_name, file_name = os.path.split(filePath)
        base_name, ext = os.path.splitext(file_name)
        prog = None
        if (ext == '.gz' or ext == '.gzip' or ext == '.bz2' or ext == '.bzip2' or ext == '.zip'):
           prog = ext[1:]
           base_name, ext = os.path.splitext(base_name)
        if len(ext) > 1:
            ext = ext[1:]
        return (dir_name, file_name, base_name, ext, prog)


    def readNameByReadFile(self, readFilePath):
        """
        get read name by read file path.
        """
        dir_name, file_name, read_name, ext, prog = self.splitFilePath(readFilePath)
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
        dir_name, file_name, basename, ext, prog = self.splitFilePath(readFile)
        if secondReadType:
            ext = ".%s" % secondReadType
        if prog is not None:
            fpath = "%s/%s2%s.%s" % (dir_name, basename[0:-1], ext, prog)
        else:
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


    def mcContextType(self, genomeSeq, cBasePos, strand='+'):
        """
        get mC context type (CG, CHG, CHH) by genome sequence and C base position.
        if no mC context found, return None.
        """

        try:
            if strand == '+':
                if genomeSeq[cBasePos + 1] == "G":
                    return "CG"
                else:
                    if genomeSeq[cBasePos + 2] == "G":
                        return "CHG"
                    else:
                        return "CHH"
            elif strand == '-':
                if genomeSeq[cBasePos - 1] == "C":
                    return "CG"
                else:
                    if genomeSeq[cBasePos - 2] == "C":
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


    def bzip2File(self, filePath, wait = True, log = False):
        """
        bzip2 file. If wait argument is False, without waiting for bzip2 process to be
        completed, this function returns immediately.
        """

        if log:
            logging.info("bzip2 start: %s" % filePath)

        dirpath, fname = os.path.split(filePath)
        cmd = "bzip2 %s" % fname
        p = subprocess.Popen(cmd, shell = True, cwd = dirpath)
        if wait:
            p.wait()

        if log:
            logging.info("bzip2 done: %s" % filePath)


    def isGzipFile(self, filePath):
        return filePath[-3:] == ".gz" or filePath[-5:] == ".gzip"


    def isBzip2File(self, filePath):
        return filePath[-4:] == ".bz2" or filePath[-6:] == ".bzip2"


    def isZipFile(self, filePath):
        return filePath[-4:] == ".zip"


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


    def readRefGenome(self, refGenome, refGenomeBuf, refGenomeChr):
        """
        read the reference genome fasta file.
        """

        logging.info("BsfCallBase::readRefGenome: %s" % refGenome)
        chr = None
        buf = []
        fin = open(refGenome, 'r')
        for line in fin:
            if line[0] == '>':
                chr = self.chrnoFromFastaDescription(line)
                logging.info("BsfCallBase::readRefGenome: chr=%s" % chr)
                if len(buf) > 0:
                    refGenomeBuf[refGenomeChr[-1]]=''.join(buf)
                    del buf[:]
                refGenomeChr.append(chr)
            elif chr != None:
                buf.append(line.strip().upper())
            else:
                logging.fatal("BsfCallBase::readRefGenome: the specified reference genome file \"%s\" is malformed." % refGenome)
        fin.close()
        refGenomeBuf[refGenomeChr[-1]]=''.join(buf)
        logging.info("BsfCallBase::readRefGenome: done.")
        return


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
                    logging.info("McDetector::getAllMappingResultFiles: %s" % filename)
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

        # self.numReadsPerFile = self.sizeForSplitRead(self.opts["split_read_size"])

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
                self.prepareForReads()
                self.mappingResultDirs = self.processMapping()

            logging.debug("BsfCall:execute: mapping result directories are: %s" % ','.join(self.mappingResultDirs))
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
            self.runLast(read_attr)
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
        # log_level = logging.DEBUG

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
            # logging.info("  Read BAM file: %s" % ("Yes" if self.opts["read_bam"] else "No"))
            # logging.info("  Read SAM file: %s" % ("Yes" if self.opts["read_sam"] else "No"))
        else:
            logging.info("Reference genome: %s" % self.refGenome)
            logging.info("Read files: %s" % self.readFilePaths)
            logging.info("Working directory: %s" % self.dataDir)
            logging.info("Options:")
            logging.info("  Threshold of the alignment score at filtering: %d" % self.opts["aln_score_thres"])
            # logging.info("  Paired-end direction: %s" % self.opts["pe_direction"])
            # logging.info("  Options for LAST: %s" % self.opts["last_opts"])

        logging.info("  Threshold of read coverate: %d" % self.opts["coverage"])
        logging.info("  Threshold of mC ratio: %s" % str(self.opts["lower_bound"]))
        logging.info("  Threshold of the mismap probability at filtering: %s" % str(self.opts["aln_mismap_prob_thres"]))
        logging.info("  Working directory: %s" % self.dataDir)
        logging.info("  Local directory: %s" % self.opts["local_dir"])
        logging.info("  Output file: %s" % (self.opts["output"] if self.opts["output"] else "(stdout)"))
        # logging.info("  Use cluster: %s" % ("Yes" if self.opts["use_cluster"] else "No"))
        # logging.info("  Queue name: %s" % self.opts["queue_list"])
        logging.info("  Number of threads: %d" % self.opts["num_threads"])
        # logging.info("  Split read size: %s" % self.opts["split_read_size"])


    def prepareForReads(self):
        """
        create directories to store split reads and result files.
        """

        for read_no, read_path in enumerate(self.readFilePaths):
            readNo = read_no + 1
            logging.info("Preparations for a read file start: %d: %s" % (readNo, read_path))

            data_dir = "%s/%d" % (self.dataDir, readNo)
            read = {"base_dir": data_dir, "path": read_path, "reads_dir": data_dir + "/reads", "results_dir": data_dir + "/results"}

            if not os.path.exists(data_dir):
                os.makedirs(data_dir)
                os.makedirs(read["reads_dir"])
                os.makedirs(read["results_dir"])

            pe_no = self.pairedEndReadNumbers()[0]
            for readpath in read_path.split(","):
                dir_name, file_name, base_name, ext, prog = self.splitFilePath(readpath)
                file_type = self.checkReadFileType(readpath)
                read[pe_no] = {"name": base_name, "fname": file_name, "type": file_type, "path": readpath}
                pe_no += 1

            is_paired_end = self.isPairedEnd(read)
            logging.info("Paired-end: %s" % is_paired_end)
            logging.info("  Forward: %s" % read[1]["path"])
            if is_paired_end:
                logging.info("  Reverse: %s" % read[2]["path"])

            logging.info("Preparations for a read file done")
            self.reads.append(read)
        return


    def runLast(self, readAttr):
        """
        run LAST programs to map reads and filtering.
        """

        is_paired_end = self.isPairedEnd(readAttr)
        filter_option = self.filterOpts(self.opts["aln_mismap_prob_thres"], self.opts["aln_score_thres"], is_paired_end)

        if is_paired_end:
            logging.info('BsfCall::runLast: PairedEnd')
            last_exec = LastExecutorPairedEnd(self.refGenome, self.dataDir, readAttr["reads_dir"], readAttr["results_dir"], self.opts["num_threads"])
        else:
            logging.info('BsfCall::runLast: Not PairedEnd')
            last_exec = LastExecutorSingle(self.refGenome, self.dataDir, readAttr["reads_dir"], readAttr["results_dir"])

        # last_exec.execute(readAttr, self.opts["num_threads"], self.lastalOpts(self.opts["last_opts"]), self.mergeOpts(), filter_option)
        # last_exec.execute(readAttr, 1, self.lastalOpts(self.opts["last_opts"]), self.mergeOpts(), filter_option)
        last_exec.execute(readAttr, 1, "", self.mergeOpts(), filter_option)


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

        # if self.opts["read_bam"]:
        #     argv.append("BAM")
        # elif self.opts["read_sam"]:
        #     argv.append("SAM")
        # else:
        #     argv.append("MAF")
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
            cmd = "qsub -o %s -e %s -q %s -pe openmpi %d -cwd -b y %s %s" % (out_file, err_file, self.opts["queue_list"], self.opts['num_threads'], remote_cmd, cmdArgs)
        else:
            cmd = "qsub -o %s -e %s -pe openmpi %d -cwd -b y %s %s" % (out_file, err_file, self.opts['num_threads'], remote_cmd, cmdArgs)

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

    def __init__(self, refGenome, baseDir = ".", readsDir = None, resultsDir = None, numThreads = 1):
        """
        constructor of LastExecutor
        """

        self.refGenome = refGenome
        self.baseDir = baseDir
        self.readsDir = readsDir
        self.resultsDir = resultsDir
        self.queue = Queue.Queue()
        self.lock = threading.Lock()
        self.numThreads = numThreads


    def execute(self, readAttr, numThreads = 1, lastalOpts = "", mergeOpts = "", filterOpts = ""):
        """
        enqueue all splited read files to the queue.
        create and start threads to execute read mapping and filtering process.
        wait for all threads to finish.
        """

        self.enqueue(readAttr)

        if self.queue.qsize()==0:
            logging.fatal("LastExecutor::execute: Error: queue size=%d" % self.queue.qsize())
            sys.exit(1)

        logging.info("Queued %d files." % self.queue.qsize())

        threads = []
        for i in range(numThreads):
            t = threading.Thread(target=self.worker, args=(readAttr, lastalOpts, mergeOpts, filterOpts))
            t.daemon = True
            threads.append(t)
            t.start()

        i = 1
        for thread in threads:
            logging.info("Waiting for %d th thread." % i)
            thread.join()
            logging.info("Joined %d th thread." % i)
            i+=1


    def worker(self, readAttr, lastalOpts, mergeOpts, filterOpts):
        """
        thread worker to execute LAST.
        dequeue read file path from the queue and execute read mapping filtering
        process.
        """
        if self.queue.empty():
            logging.info("LastExecutor::worker: Queue is empty.")
        while not self.queue.empty():
            fpath = self.queue.get_nowait()
            self.runLast(fpath, readAttr, lastalOpts, mergeOpts, filterOpts)
        return

        
    def runLast(self, readFile, readAttr, lastalOpts, mergeOpts, filterOpts, rmInFiles = True):
        """
        execute LAST programs to map read and filtering.
        """

        cmd = self.batchCmd(readFile, readAttr, lastalOpts, mergeOpts, filterOpts, rmInFiles)
        logging.info("LastExecutor::runLast: command=%s" % cmd)
        p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        # out, error = p.communicate()
        p.wait()

        error_msg = p.communicate()[1]
        if len(error_msg) > 0:
            logging.fatal(error_msg)


    def enqueue(self, readAttr):
        """
        enqueue all splitted read files to the queue.
        """

        # logging.info("%s/*_1.%s.bz2" % (self.readsDir, readAttr[1]["type"]))
        # for read_file in glob.glob("%s/*_1.%s" % (self.readsDir, readAttr[1]["type"])):
        #     self.queue.put(read_file)
        self.queue.put(readAttr[1]["path"])


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
        # return "lastdb -w2 -u bisulfite_%s.seed %s.%s %s" % (direction, self.refGenome, direction, self.refGenome)
        return "lastdb -uBIS%s %s.%s %s" % (direction.upper(), self.refGenome, direction, self.refGenome)


    def lastalBsCmd(self, readFile, opts = ""):
        """
        get lastal command to map read.
        """
        # s_opt = self.lastalSopt(direction)
        read_name = self.readNameByReadFile(readFile)
        return "TMPDIR=%s last-bisulfite.sh %s.%s %s.%s %s > %s" % (self.readsDir, self.refGenome, 'f', self.refGenome, 'r', readFile, self.mappingResultFilePath(read_name))

    def lastalBsPairCmd(self, readFile1, readFile2, opts = ""):
        """
        get lastal command to map read.
        """
        # s_opt = self.lastalSopt(direction)
        read_name = self.readNameByReadFile(readFile)
        return "PARALLEL=-j%d TMPDIR=%s last-bisulfite-pair.sh %s.%s %s.%s %s %s > %s" % (self.numThreads, self.readsDir, self.refGenome, 'f', self.refGenome, 'r', readFile1, readFile2, self.mappingResultFilePath(read_name))


    # def mergeCmd(self, forwardFile, reverseFile, outputFile, opts = "", rmInFiles = True):
    def mergeCmd(self, inputFile, outputFile, opts = "", rmInFiles = True):
        """
        get command to merge lastal output.
        """
        cmd = "last-merge-batches %s > %s" % (inputFile, outputFile)
        if rmInFiles:
            cmd += "; rm %s" % inputFile
        return cmd


    # def mappingAndMergeCmd(self, readFile, lastalOpts = "", mergeOpts = "", rmInFiles = True):
    def mappingPairCmd(self, readFile1, readFile2, lastalOpts = "", mergeOpts = "", rmInFiles = True):
        """
        get read mapping and filtering command.
        """
        read_name = self.readNameByReadFile(readFile)
        n, ext = os.path.splitext(readFile)
        if ext == ".gz":
            n, ext = os.path.splitext(n)
        lastal_qopt = self.lastalQopt(ext[1:])
        lastal_opt = "%s %s" % (lastalOpts, lastal_qopt)
        mapping_file = self.mappingResultFilePath(read_name)

        return self.lastalBsPairCmd(readFile1, readFile2, lastal_opt)


    # def mappingResultFilePath(self, readName, direction):
    def mappingResultFilePath(self, readName):
        """
        get read mapping result file path.
        """

        return "%s/%s.maf" % (self.resultsDir, readName)


    def mergeResultFilePath(self, readName):
        """
        get merge result file path.
        """

        return "%s/%s.merge.maf" % (self.resultsDir, readName)


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

        dir_path, file_name, base_name, ext, prog = self.splitFilePath(readFile)
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


    def batchCmd(self, readFile1, readAttr, lastalOpts = "", mergeOpts = "", filterOpts = "", rmInFiles = True):
        """
        get batch command to map read and filtering.
        """
        read_name = self.readNameByReadFile(readFile1)
        out_file = self.filterResultFilePath(read_name[0:-2])
        return "TMPDIR=%s last-bisulfite.sh %s.f %s.r %s > %s/%s.maf" % (self.readsDir, self.refGenome, self.refGenome, readFile1, self.resultsDir, read_name)
        # cmds = []
        # cmds.append(self.mappingAndMergeCmd(readFile, lastalOpts, mergeOpts, rmInFiles))
        # cmds.append(self.mappingPairCmd(readFile1, readFile2, lastalOpts, mergeOpts, rmInFiles))
        # cmds.append(self.filterCmd(self.mergeResultFilePath(read_name), out_file, filterOpts, rmInFiles))
        # cmds.append(self.filterCmd(self.mappingResultFilePath(read_name), out_file, filterOpts, rmInFiles))
        # cmds.append("bzip2 %s" % out_file)
        # return "; ".join(cmds)


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

    def __init__(self, refGenome, baseDir, readsDir, resultsDir, numThreads):
        """
        constructor of LastExecutorPairedEnd
        """
        LastExecutor.__init__(self, refGenome, baseDir, readsDir, resultsDir)


    def batchCmd(self, readFile1, readAttr, lastalOpts = "", mergeOpts = "", filterOpts = "", rmInFiles = True):
        """
        get batch command to map read and filtering.
        """
        # readFile2 = self.secondReadFilePathByFirstReadFilePath(readFile1, readAttr[2]["type"])
        # read_name1 = self.readNameByReadFile(readFile1)
        # read_name2 = self.readNameByReadFile(read_files[1])
        # merge_result_file = "%s %s" % (self.mergeResultFilePath(read_name1), self.mergeResultFilePath(read_name2))
        # mapping_result_file = "%s %s" % (self.mappingResultFilePath(read_name1), self.mappingResultFilePath(read_name2))
        # out_file = self.filterResultFilePath(read_name1[0:-2])
        # return "TMPDIR=%s last-bisulfite-paired.sh %s.f %s.r %s %s > %s" % (self.baseDir, self.refGenome, self.refGenome, readFile1, readFile2, out_file)
        read_name1 = self.readNameByReadFile(readAttr[1]["path"])
        read_name2 = self.readNameByReadFile(readAttr[2]["path"])
        return "PARALLEL=-j%d TMPDIR=%s last-bisulfite-paired.sh %s.f %s.r %s %s > %s/%s,%s.maf" % (self.numThreads, self.readsDir, self.refGenome, self.refGenome, readAttr[1]["path"], readAttr[2]["path"], self.resultsDir, read_name1, read_name2)


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
        self.refGenomeBuf = {}
        self.refGenomeChr = []

        self.mappingResultDirs = resultDirs
        self.mcContextDir = mcContextDir
        self.lowerBound = options["lower_bound"]
        self.coverageThreshold = options["coverage"]
        self.onlyMcDetection = options["only_mcdetection"]

        self.opts = options

        self.mappingResultFiles = []

        self.mismapThreshold = options["aln_mismap_prob_thres"]
        self.readBam = False
        self.readSam = False

        if self.onlyMcDetection:
            if "mapping_result_files" in options:
                self.mappingResultFiles = options["mapping_result_files"]
            else:
                self.mappingResultFiles = self.getAllMappingResultFiles(resultDirs)
        else:
            self.mappingResultFiles = self.getAllMappingResultFiles(resultDirs)

        if len(self.mappingResultFiles)==0:
            logging.fatal("McDetector::__init__: error: no mapping result file found.")
            sys.exit(1)
        
        if options["local_dir"]:
            self.localDir = options["local_dir"]
        else:
            self.localDir = mcContextDir

        self.readRefGenome(self.refGenome, self.refGenomeBuf, self.refGenomeChr)

        logging.debug("McDetector::__init__: mappingResultDirs=%s" % ','.join(self.mappingResultDirs))


    def execute(self, outputFile, numWorkers = 1):
        """
        execute mC detection process and output result.
        """

        self.process()
        self.output(outputFile)


    def process(self):
        """
        execute mC detection process.
        """

        logging.info("mC detection process start")

        if len(self.mappingResultFiles)==0:
            logging.fatal("McDetector::process: error: no mapping result found.")
            sys.exit(1)

        for mapping_result_file in self.mappingResultFiles:
            logging.info("Parsing mapping result file: %s" % mapping_result_file)
            if self.onlyMcDetection:
                if not (self.readBam or self.readSam):
                    if self.isGzipFile(mapping_result_file) or self.isMafFile(mapping_result_file):
                        self.processMafFile(mapping_result_file)
                else:
                    # BAM or SAM
                    if self.readBam and self.isBamFile(mapping_result_file):
                        self.processSamFile(mapping_result_file)
                    elif self.readSam and self.isSamFile(mapping_result_file):
                        self.processSamFile(mapping_result_file)
            else:
                self.processMafFile(mapping_result_file)

        logging.info("Parsing mapping result file done")

        if self.mcContextDir != self.localDir:
            copy(self.mcContextLocalFilePath(self.targetChr), self.mcContextDir)


    def output(self, outputFile):
        """
        output mC detection result.
        """

        logging.info("McDetector::output: outputFile=%s" % outputFile)
        popen_args = ['sort', '-k1']
        list = glob.glob("%s/*.bsf" % self.localDir)
        if len(list)==0:
            logging.fatal("McDetect::output: no *._bsf_ files found in %s" % self.localDir)
            sys.exit(1)
        for bsf_file in list:
            if os.path.getsize(bsf_file) > 0:
                popen_args.append(bsf_file)
                logging.info("McDetector::output: added bsf_file=\"%s\"" % bsf_file)
        logging.debug("McDetector::output: popen_args=%s" % ' '.join(popen_args))
        fout = open(outputFile, 'w')
        if len(popen_args) > 2:
            pipe = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
            block = []
            for line in pipe.stdout:
                (chr, ctx_pos, strand, mc_ctx, base, conf) = line.strip().split("\t")
                ctx_pos = int(ctx_pos)
                conf = float(conf)
                if len(block)==0 or block[-1][0]==chr:
                    block.append([chr, ctx_pos, strand, mc_ctx, base, conf])
                else:
                    self.outputOneChrBlock(fout, block)
                    del block[:]
                    block.append([chr, ctx_pos, strand, mc_ctx, base, conf])
            if len(block) > 0:
                self.outputOneChrBlock(fout, block)
            pipe.stdout.close()
        else:
            logging.info("McDetector::output: no result files found.")
        fout.close()


    def outputOneChrBlock(self, fout, block):
        """
        output mC detection result for one chromosome.
        """

        chr = block[0][0]
        logging.info("McDetector::outputOneChrBlock: chr=%s" % chr)
        mc_contexts = {}
        for b in sorted(block, key=lambda block: block[1]):
            try:
                ctx_pos, strand, mc_ctx, base, conf = b[1:]
                if not ctx_pos in mc_contexts:
                    mc_contexts[ctx_pos] = {}
                if not strand in mc_contexts[ctx_pos]:
                    mc_contexts[ctx_pos][strand] = {}
                if not mc_ctx in mc_contexts[ctx_pos][strand]:
                    mc_contexts[ctx_pos][strand][mc_ctx] = []
                mc_contexts[ctx_pos][strand][mc_ctx].append([base, conf])
            except ValueError, e:
                logging.warning("McDetect::outputOneChrBlock: value error: %s: %s -> %s" % (fpath, line.strip(), e.args[0]))

        num_entry = 0
        if len(mc_contexts.keys()) > 0:
            for pos in sorted(mc_contexts.keys()):
                for strand in mc_contexts[pos].keys():
                    for mc_ctx in mc_contexts[pos][strand].keys():
                        coverage, mc_ratio = self.calcCoverage(mc_contexts[pos][strand][mc_ctx], strand)
                        if coverage >= self.coverageThreshold and mc_ratio >= self.lowerBound:
                            fout.write("%s\t%d\t%s\t%s\t%g\t%d\n" % (chr, pos, strand, mc_ctx, mc_ratio, coverage))
                            num_entry += 1
                        else:
                            logging.info("McDetect::outputOneChrBlock: rejected: chr=%s pos=%d strand=%s coverage=%g mc_ratio=%g" % (chr, pos, strand, coverage, mc_ratio))
            # self.bzip2File(fpath, False)
        logging.info("McDetector::outputOneChrBlock: number of entries %s" % num_entry)


    def processMafFile(self, resultFile):
        """
        read mapping result file, and output mC context data file.
        """

        file_name = self.splitFilePath(resultFile)[1]
        outputFile = "%s/%s.bsf" % (self.localDir, file_name)
        try:
            fin = open(resultFile, 'r')
            fout = open(outputFile, 'w')
            block = []
            for line in fin:
                if line[0] == '#' or line[0] == 'p' or line[0] == "\n":
                    continue
                if line[0] == 'a' or line[0] == 's' or line[0] == 'q':
                    block.append(line.strip())
                if len(block)==4:
                    if block[0][0]=='a' and block[1][0]=='s' and block[2][0]=='s' and block[3][0]=='q':
                        mismap_prob = float(block[0].split('=')[2])
                        if mismap_prob <= self.mismapThreshold:
                            b1 = block[1].split()
                            b2 = block[2].split()
                            chr = b1[1]
                            ref_seq = b1[-1].upper()
                            ref_start = int(b1[2]) # 0-based position
                            ref_len = int(b1[3])
                            ref_size = int(b1[5])
                            read_seq = b2[-1]
                            read_len = int(b2[3])
                            read_qual = block[3].split()[-1]
                            strand = self.findStrandFromAlignment(ref_seq, read_seq)
                            # strand = block[2].split()[4]
                            if strand == '+' or strand == '-':
                                lines = self.extractMcContextsByOneRead(chr, strand, mismap_prob, ref_seq, ref_start, read_seq, read_qual, read_len)
                                for line in lines:
                                    fout.write(line)
                                logging.debug("processMafFile: a maf block(%s) is successfully captured." % strand)
                            elif strand == '+-':
                                for st in ('+', '-'):
                                    lines = self.extractMcContextsByOneRead(chr, st, mismap_prob, ref_seq, ref_start, read_seq, read_qual, read_len)
                                    for line in lines:
                                        fout.write(line)
                                logging.debug("processMafFile: a maf block(%s) is successfully captured." % strand)
                            else:
                                logging.debug("processMafFile: a maf block does not show strand-specific info.")
                        else:
                            logging.debug("processMafFile: alignment \"%s\" has greater mismap prob. than the threshold." % block[0])
                        del block[:]
                    else:
                        logging.fatal("processMafFile: error: unexpected malformed maf block is found.")
                        logging.fatal("block 1: \"%s\"\n" % block[0])
                        logging.fatal("block 2: \"%s\"\n" % block[1])
                        logging.fatal("block 3: \"%s\"\n" % block[2])
                        logging.fatal("block 4: \"%s\"\n" % block[3])
                        sys.exit(1)
            fin.close()
            fout.close()

            if len(block) > 0:
                logging.fatal("McDetect::processMafFile: error: possible malformed MAF file.")
                for b in block:
                    logging.fatal(b)
                sys.exit(1)
        except IOError:
            logging.fatal("McDetect::processMafFile: error: unable to read a MAF file \"%s\"" % resultFile)
            sys.exit(1)
        return
     

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
                self.extractMcContextsByOneRead(chr, samaln.strand, 0, samaln.alnRefSeq.upper(), samaln.alnRefStart, self.targetSeqLen, read_seq, read_qual)

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

    def extractMcContextsByOneRead(self, chr, strand, mismapProb, refSeq, refStart, readSeq, readQual, readLen):
        """
        extract mC context by one read.
        """

        lines = []
        if strand == '+':
            real_pos = refStart
            for i in range(readLen):
                if readSeq[i] == '-':
                    continue
                # if refSeq[i] == 'C' and (readSeq[i] == 'C' or readSeq[i] == 'T'):
                if refSeq[i] == 'C':
                    baseConf = self.qualityCharToErrorProb(readQual[i])
                    ctx_type = self.mcContextType(self.refGenomeBuf[chr], real_pos, strand)
                    line = "%s\t%d\t%s\t%s\t%s\t%g\n" % (chr, real_pos, strand, ctx_type, readSeq[i], (1-mismapProb) * (1-baseConf))
                    # logging.debug("extractMcContextsByOneRead: line=%s" % line)
                    lines.append(line)
                real_pos += 1
        if strand == '-':
            real_pos = refStart
            for i in range(readLen):
                if readSeq[i] == '-':
                    continue
                # if refSeq[i] == 'G' and (readSeq[i] == 'G' or readSeq[i] == 'A'):
                if refSeq[i] == 'G':
                    baseConf = self.qualityCharToErrorProb(readQual[i])
                    ctx_type = self.mcContextType(self.refGenomeBuf[chr], real_pos, strand)
                    read_base = self.complementaryBase(readSeq[i])
                    line = "%s\t%d\t%s\t%s\t%s\t%g\n" % (chr, real_pos, strand, ctx_type, readSeq[i], (1-mismapProb) * (1-baseConf))
                    # logging.debug("extractMcContextsByOneRead: line=%s" % line)
                    lines.append(line)
                real_pos += 1
        return lines


    def complementaryBase(self, base):
        """
        compute a complemetary base.
        """

        if base == 'A': return 'T'
        if base == 'C': return 'G'
        if base == 'G': return 'C'
        if base == 'T': return 'A'
        if base == 'N': return 'N'
        if base == 'a': return 't'
        if base == 'c': return 'g'
        if base == 'g': return 'c'
        if base == 't': return 'a'
        if base == 'n': return 'n'
        sys.exit(1)


    def findStrandFromAlignment(self, ref_seq, read_seq):
        plus_sup = 0
        minus_sup = 0
        other_mismatch = 0
        base_size = 0
        for i in range(len(read_seq)):
            if ref_seq[i]=='C' and read_seq[i]=='T':
                plus_sup += 1
                base_size += 1
            elif ref_seq[i]=='G' and read_seq[i]=='A':
                minus_sup += 1
                base_size += 1
            elif (ref_seq[i]!='-' and read_seq[i]!='-') and ref_seq[i]!=read_seq[i]:
                other_mismatch += 1
                base_size += 1
        if base_size==0:
            return '.'
        mismatch_rate = float(other_mismatch)/float(base_size)
        # if plus_sup > minus_sup and plus_sup > other_mismatch:
        # if plus_sup > minus_sup and mismatch_rate < 0.05:
        if plus_sup > minus_sup:
            return '+'
        # if minus_sup > plus_sup and minus_sup > other_mismatch:
        # if minus_sup > plus_sup and mismatch_rate < 0.05:
        if minus_sup > plus_sup:
            return '-'
        if plus_sup == minus_sup:
            return '+-'
        return '.'


    def __extractMcContextsByOneRead(self, chr, strand, mismapProb, refSeq, refStart, refSeqLen, readSeq, readQual):
        """
        extract mC context by one read.
        """

        logging.debug("extractMcContextsByOneRead(%s, %s, %g, %s, %d, %d, %s, %s)" % (chr, strand, mismapProb, refSeq, refStart, refSeqLen, readSeq, readQual))

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
                if ctx_type == None:
                    continue
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

    def qualityCharToErrorProb(self, qualityChar):
        """
        convert FASTQ (Sanger) quality character to error probability.
        """
        return 10**((ord('!')-ord(qualityChar))*0.1)

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


    def calcCoverage(self, seqAry, strand):
        """
        count the number of C and T (or G and A), calculate mC ratio.
        """

        num_c = 0.0
        num_t = 0.0
        num_all = 0
        (C, T) = ('C', 'T')
        if strand == '-':
            (C, T) = ('G', 'A')
        for i in seqAry:
            num_all += 1
            if i[0] == C:
                num_c += i[1] 
            elif i[0] == T:
                num_t += i[1]

        num_ct = num_c + num_t
        if num_all == 0 or num_ct == 0:
            return (0, 0)
        # return (num_all, num_c/num_all)
        return (num_all, num_c/num_ct)


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


    def mcContextLocalFilePath(self, chrNo):
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

