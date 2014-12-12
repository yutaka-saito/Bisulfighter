Bisulfighter: a pipeline for accurate detection of methylated cytosines and differentially methylated regions 

Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
by National Institute of Advanced Industrial Science and Technology (AIST)
is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
http://creativecommons.org/licenses/by-nc-sa/3.0/

Toutai Mituyama <mituyama-toutai at aist dot go dot jp>
Computational Biology Research Center
National Institute for Advanced Industrial Science and Technology (AIST)

--------------------------------------------
I. PREREQUISITS

1) LAST alignment program
    Please get it from the original distribution at http://last.cbrc.jp/

2) Python (2.4.x)
    Bisulfighter scripts are written in Python.

3) R (optional)
    bsf-diff uses R.

4) Boost C++ library
    ComMet requires Boost C++ library.
    http://www.boost.org/


II. PACKAGE COMPONENTS

1) bsf-call
    One of the major components of Bisulfighter, which is a python script for mC-call pipeline perfoms:
    a) short read mapping with LAST
    b) mC detection and mC rate estimation

2) ComMet
    A HMM-based differentially methylated region (DMR) identifier.

3) bsf-diff
    An optional component for DMR identification using statistical test and signal smoothing.
    This program is not officially released. You can use this one but we are not going to provide details about this.

4) script
    Scripts to generate simulated reads for differentially methylated
cytosines and differentially methylated reagions.

5) demo
    Scripts to perform a simple demonstration of bsf-call and ComMet.


III. INSTALLATION

1) bsf-call
    You can place it anywhere you want but bsf-call comsumes great amount of
disk sapce for its working directory. So it is better to place it on a
high-performance large disk volume.
    bsf-call comsumes large physical memory space while it perfoms LAST
alignments because its database resides on-memory entirely.
    For example, 30GB or more is recommended for hg19 human genome.

2) ComMet
     You can place its executable "ComMet" anywhere you want.

EOT

