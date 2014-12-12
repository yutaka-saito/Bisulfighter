Bisulfighter: a pipeline for accurate detection of methylated cytosines and differentially methylated regions 

Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
by National Institute of Advanced Industrial Science and Technology (AIST)
is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0
Unported License. http://creativecommons.org/licenses/by-nc-sa/3.0/

--------------------------------------------
I. PREREQUISITS

1) LAST alignment program
    Please get it from the original distribution at http://last.cbrc.jp/

2) Python (2.4.x)
    Bisulfighter scripts are written in Python.

3) Boost C++ library 
    ComMet uses Boost C++ library
    http://www.boost.org/

II. PACKAGE COMPONENTS

1) bsf-call
    One of the major components of Bisulfighter, which is a python script
    for mC-call pipeline perfoms:
    a) short read mapping with LAST
    b) mC detection and mC rate estimation

2) ComMet
    A HMM-based differentially methylated region (DMR) identifier.

III. INSTALLATION
    cd ComMet
    make


IV. DEMO
    cd demo
    ./demo.sh 2>&1 | tee demo.log


