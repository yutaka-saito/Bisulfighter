#!/usr/bin/env python

import sys
import os
import logging
import subprocess
import bsfcall

__version__ = "0.1"

refgenome = sys.argv[1]
mapping_dirs = sys.argv[2].split(',')
work_dir = sys.argv[3]
chr_no = sys.argv[4]
if sys.argv[5] == "True":
    only_mcdetection = True
else:
    only_mcdetection = False

try:
    lower_bound = float(sys.argv[6])
except Exception:
    lower_bound = 0.01

try:
    coverage_threshold = int(sys.argv[7])
except Exception:
    coverage_threshold = 5

try:
    mismap = float(sys.argv[8])
except Exception:
    mismap = 0.01

try:
    target_format = sys.argv[9]
except:
    target_format = "MAF"

try:
    local_dir = (sys.argv[10])
except Exception:
    local_dir = None

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

if local_dir and not os.path.exists(local_dir):
    os.makedirs(local_dir)


log_level = logging.INFO
log_file = "%s/mc-detector-%s.log" % (work_dir, chr_no)
file_logger = logging.FileHandler(filename=log_file)

file_logger.setLevel(log_level)
file_logger.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))

logging.getLogger().addHandler(file_logger)
logging.getLogger().setLevel(log_level)

hostname = subprocess.check_output(['hostname'])

logging.info("%s start." % sys.argv[0])
logging.info("Host: %s" % hostname.strip())
logging.info("Arguments:")
logging.info("  Reference genome: %s" % refgenome)
logging.info("  Mapping result directory: %s" % mapping_dirs)
logging.info("  Only mC detection: %s" % str(only_mcdetection))
logging.info("  Working directory: %s" % work_dir)
logging.info("  Chromosome: %s" % chr_no)
logging.info("  Threshold of read coverate: %d" % coverage_threshold)
logging.info("  Threshold of mC ratio: %s" % str(lower_bound))
logging.info("  Threshold of the mismap probability at filtering: %s" % str(mismap))
logging.info("  Target mapping file format: %s" % target_format)
logging.info("  Local directory: %s" % local_dir)

options = {}
options["only_mcdetection"] = only_mcdetection
options["lower_bound"] = lower_bound
options["coverage"] = coverage_threshold
options["aln_mismap_prob_thres"] = mismap
options["local_dir"] = local_dir

options["read_bam"] = False
options["read_sam"] = False
options["bam2sam_dir"] = None

if target_format == "SAM":
    options["read_sam"] = True
elif target_format[0:3] == "BAM":
    options["read_bam"] = True
    options["bam2sam_dir"] = target_format[4:]

mc_detector = bsfcall.McDetector(refgenome, mapping_dirs, work_dir, options)
mc_detector.processOneChr(chr_no)

logging.info("%s done." % sys.argv[0])

sys.exit(0)
