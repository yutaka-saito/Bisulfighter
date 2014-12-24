#!/bin/sh
# 
# Bisulfighter Demonstration Script
# 
# Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
# Usage: ./demo.sh 2>&1 | tee demo.log

invoke_command () {
  echo "invoking: $1"
  eval $1
}

remove_if () {
  for i in $1 ; do
    if [ -d $i -o -f $i ] ; then
      invoke_command "\rm -rf $i"
    fi
  done
}

#
# The following lines work unexpectedly in sun grid engine environment.
# Please be sure to run this script in demo/ directory as your currently working directory
#
# TOPDIR=`dirname $0`
# cd $TOPDIR

typeset -i PYTHON_VERSION=`/usr/bin/env python -c 'import sys; print "%d" % sys.hexversion'`
if [ $PYTHON_VERSION -lt 34013184 ];
then
  echo "Error: bsf-call requires Python-2.7 or later."
  exit 1;
fi
BSFCALL_SCRIPT=../bsf-call/bsf-call
if [ ! -f $BSFCALL_SCRIPT ];
then
  echo "Error: ${BSFCALL_SCRIPT} is not found."
  exit 1;
else
  ln -fs $BSFCALL_SCRIPT .
  BSFCALL_SCRIPT=./bsf-call
fi;
BSFCALL_PY=../bsf-call/bsfcall.py
if [ ! -f $BSFCALL_PY ];
then
  echo "Error: $BSFCALL_PY is not found."
  exit 1;
else
  ln -fs $BSFCALL_PY .
fi;
MAPPING_P_SH=../bsf-call/mapping-p.sh
if [ ! -f $MAPPING_P_SH ];
then
  echo "Error: $MAPPING_P_SH is not found."
  exit 1;
else
  ln -fs $MAPPING_P_SH .
fi;
MAPPING_S_SH=../bsf-call/mapping-s.sh
if [ ! -f $MAPPING_S_SH ];
then
  echo "Error: $MAPPING_S_SH is not found."
  exit 1;
else
  ln -fs $MAPPING_S_SH .
fi;
MC_DETECTOR_PY=../bsf-call/mc-detector.py
if [ ! -f $MC_DETECTOR_PY ];
then
  echo "Error: $MC_DETECTOR_PY is not found."
  exit 1;
else
  ln -fs $MC_DETECTOR_PY .
fi;


BSFCALL_CMD="${BSFCALL_SCRIPT} -c 0 --last=-d108,-e120"
BSFCALL_README=../bsf-call/bsf-call.txt
COMMET_DIR=../ComMet
COMMET_EXEC=${COMMET_DIR}/src/ComMet
if [ ! -f $COMMET_EXEC ];
then
  echo "${COMMET_EXEC} is not found."
  if [ -d $COMMET_DIR ];
  then
    echo "Build ComMet executable."
    ( cd $COMMET_DIR; make )
    if [ -x $COMMET_DIR ];
    then
      echo "done"; echo;
    else
      echo "Error: Building ComMet failed."
      exit 1;
    fi;
  fi;
fi;
COMMET_CMD="${COMMET_EXEC}"
COMMET_README=../ComMet/README
BSF2COM_SCRIPT=../ComMet/util/Bsf2ComMetIn.pl
if [ ! -f $BSF2COM_SCRIPT ];
then
  echo "Error: ${BSF2COM_SCRIPT} is not found."
  exit 1;
fi;
BSF2COM_CMD="${BSF2COM_SCRIPT}"

DATA=data
CHROM=$DATA/chrX.sub.fa
READ1=$DATA/1.fastq
READ2=$DATA/2.fastq
TMPMC1=tmp.mc1
TMPMC2=tmp.mc2
TMPDMR=tmp.dmr
COMIN=$TMPDMR/commet.in
RESULT=result
if [ ! -d $RESULT ];
then
  mkdir $RESULT;
fi;
RESMC1=$RESULT/res.mc1
RESMC2=$RESULT/res.mc2
RESDMR=$RESULT/res.dmr
RESDMC=$RESULT/res.dmc


echo "This is a demo script for bisulfighter."
echo ""
echo "Let us consider that we have two human samples from different biological conditions,"
echo "and want to perform mC calling and DMR detection on the subset of the chromosome X." 
echo "We have bisulfite-converted reads for each of two samples: $READ1 and $READ2,"
echo "and the reference chromosome: $CHROM"
echo ""

echo "(1) Remove temporal files for previous demo."
remove_if "$CHROM.*"
remove_if "$TMPMC1"
remove_if "$TMPMC2"
remove_if "$TMPDMR"
remove_if "$RESMC1"
remove_if "$RESMC2"
remove_if "$RESDMR"
remove_if "$RESDMC"
echo "done"; echo ""

echo "(2) Perform mC calling for sample1. (take a few minutes)"
invoke_command "$BSFCALL_CMD -o $RESMC1 -W $TMPMC1 $CHROM $READ1"
echo "done"; echo ""

echo "(3) Perform mC calling for sample2. (take a few minutes)"
invoke_command "$BSFCALL_CMD -o $RESMC2 -W $TMPMC2 $CHROM $READ2"
echo "done"; echo ""

echo "(4) Prepare input files for DMR detection."
invoke_command "mkdir $TMPDMR"
invoke_command "$BSF2COM_CMD --cpg --sample1 $RESMC1 --sample2 $RESMC2 > $COMIN"
echo "done"; echo ""

echo "(5) Perform DMR detection between sample1 and sample2."
invoke_command "$COMMET_CMD $COMIN $RESDMC $RESDMR"
echo "done"; echo ""

echo "Demo finished. Results are following:"
echo "mC calling for sample1, $RESMC1"
echo "mC calling for sample2, $RESMC2"
echo "DMR detection, $RESDMR"
echo "For more details of input/output formats and option tuning," 
echo "see $BSFCALL_README and $COMMET_README"

exit 0;

