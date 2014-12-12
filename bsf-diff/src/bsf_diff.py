#!/usr/bin/env python
# encoding: utf-8
'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''

import sys
import os

from datetime import datetime
from optparse import OptionParser
from gene import GeneParser
from threadpool import ThreadPool
from sampledata import SampleData
from significance_test import SignificanceTest
from rpy import r

__all__ = []
__version__ = 0.1
__date__ = '2012-12-13'
__updated__ = '2013-01-15'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

def preprocess_samples(sample_data):
    print sample_data
    sample_data.read_from_file()

def assign_dmr_id(dmr_dict):
    serial_no = 0
    for chrom in dmr_dict:
        for dmr in dmr_dict[chrom]:
            serial_no += 1
            dmr.dmr_id = serial_no

def associate_gene(gene_container, dmr_dict, file_name, assoc_dist, argv):
    
    try:
        try:
            file_obj = open(file_name, 'w')
            file_obj.write("# bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphyte Data\n")
            file_obj.write("# Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.\n")
            file_obj.write("# %s\n"%" ".join(argv))
            file_obj.write("# %s\n"%datetime.now().strftime('%m/%d/%Y %H:%M:%S (%Z)'))
            
            for chrom in gene_container.gene_dict.keys():
                if not chrom in dmr_dict:
                    continue
                chrom_dmr_list = dmr_dict[chrom]
                for gene in gene_container.get_of_chrom(chrom):
                    gene_start = None
                    gene_end = None
                    if gene.is_plus_strand:
                        if 0 < assoc_dist:
                            gene_start = gene.start
                            gene_end = gene.end + assoc_dist
                        else:
                            gene_start = gene.start + assoc_dist
                            gene_end = gene.end
                    else:
                        if 0 < assoc_dist:
                            gene_start = gene.start - assoc_dist
                            gene_end = gene.end
                        else:
                            gene_start = gene.start
                            gene_end = gene.end - assoc_dist
                    
                    hit_dmr_list = list()
                    for dmr in chrom_dmr_list:
                        is_dmr_start_in_gene = (gene_start <= dmr.start_pos and dmr.start_pos <= gene_end )
                        is_dmr_end_in_gene = (gene_start <= dmr.end_pos and dmr.end_pos <= gene_end )
                        is_gene_start_in_dmr = (dmr.start_pos <= gene_start and gene_start <= dmr.end_pos)
                        if is_dmr_start_in_gene or is_dmr_end_in_gene or is_gene_start_in_dmr:
                            hit_dmr_list.append(dmr)
                    
                    if 1 <= len(hit_dmr_list):
                        first_dmr = hit_dmr_list[0]
                        hit_dmr_id_list = list()
                        for hit_dmr in hit_dmr_list:
                            hit_dmr_id_list.append("DMR%07d"%hit_dmr.dmr_id)
                        file_obj.write("%s\t%f\t%f\t%f\t%f\t%s"%(gene.gene_id, first_dmr.mc_rate_1, \
                            first_dmr.mc_rate_2, first_dmr.mc_stddev_1, first_dmr.mc_stddev_2,",".join(hit_dmr_id_list)))
        except IOError, err:
            sys.stderr.write("Failed to output dmr file %s",file_name)
            sys.exit(1)
        
    finally:
        file_obj.close()

def smooth_data(data):
    sample_data=data[0]
    window_size=data[1]
    for rep_num in range(sample_data.get_number_of_replicates()):
        for chrom in sample_data.get_chromosome_list():
            met_manager = sample_data.get_manager_of_chrom(chrom)
            pos=[]
            m=[]
            cov=[]
            for methyl_c in met_manager:
                pos.append(methyl_c.position)
                m.append(methyl_c.get_methylrate(rep_num))
                cov.append(methyl_c.get_coverage(rep_num))
            r.warnings()
            r.library("locfit")
            r.assign("pos",pos)
            r.assign("m",m)
            r.assign("cov",cov)
            r.assign("h",window_size)
            r("posm=data.frame(pos,m)")
            r("fit=locfit(m~lp(pos,h=h),data=posm,maxk=1000000,weights=cov)")
            r("pp=preplot(fit,where='data',band='local',newdata=data.frame(pos=pos))")
            fit=r("pp")["fit"]
            list=r("unlist(pp$xev$xev)")
            for i, each in enumerate(list):
                position=int(each[0])
                methyl_c=met_manager.get_methyl_c(position)
                if methyl_c:
                    smoothedrate=None
                    if 1 <= fit[i]:
                        smoothedrate=1
                    elif fit[i] <= 0:
                        smoothedrate=0
                    else:
                        smoothedrate=fit[i]
                    methyl_c.update_methylrate(rep_num,smoothedrate)
                else:
                    sys.stderr.write("methyl_c doesn't exist at %d",position)
                    sys.exit(1)

def output_dmr(chrom_order, dmr_dict, file_name, argv):
    
    try:
        try:
            file_obj = open(file_name, 'w')
            file_obj.write("# bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphyte Data\n")
            file_obj.write("# Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.\n")
            file_obj.write("# %s\n"%" ".join(argv))
            file_obj.write("# %s\n"%datetime.now().strftime('%m/%d/%Y %H:%M:%S (%Z)'))
            for chrom in chrom_order:
                for dmr in dmr_dict[chrom]:
                    file_obj.write("%s\t%d\t%d\tDMR%07d\t%s\t%f\t%f\t%f\t%f\t%s\t%f\n"%(dmr.chrom, dmr.start_pos, dmr.end_pos + 1, dmr.dmr_id, dmr.strand, dmr.mc_rate_1, \
                                        dmr.mc_rate_2, dmr.mc_stddev_1, dmr.mc_stddev_2, str(dmr.fold_change), dmr.p_value))
        
        except IOError, err:
            sys.stderr.write("Failed to output dmr file %s",file_name)
            sys.exit(1)
    
    finally:
        file_obj.close()

def main(argv=None):
    '''Command line options.'''
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v1.0"
    program_build_date = "%s" % __updated__
 
    program_version_string = "%%prog %s (%s)" % (program_version, program_build_date)
    program_longdesc = "%prog [options] sample1_rep1[,sample1_rep2,...] sample2_rep1[,sample2_rep2,...] [sampleN_rep1,sampleN_rep2]"
    program_license = ""
    
    '''Variable Definitions'''
    sample_list=[]
 
    '''Store command line for later use'''
    saved_argv = list(sys.argv)
    
    ''' Parse Options'''
    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string, usage=program_longdesc, description=program_license)
        parser.add_option("-s", "--smooth", action="store_true", dest="smooth", help="smooth methylation rate data before DMR detection [default: %default]", default=False, metavar="<True|False>")
        parser.add_option("--smooth-window", type="int", dest="smooth_window", help="width of smoothing window [default: %default]", default=100, metavar="<Integer>")
        parser.add_option("-t", "--target", type="string", dest="target", help="target region in BED4 format. if this option is on, only target region is checked to be DMR or not", metavar="<BED file name>")
        parser.add_option("--gtf", type="string", dest="gtf", help="genes in GTF (Gene Transfer Format). tell bsf_diff to check DMRs are associated with these genes", metavar="<GTF file name>")
        parser.add_option("--assoc-dist", type="int", dest="assoc_dist", help="max DMR association distance for upstream and/or downstream of gene. exp) +3000", default=+3000, metavar="<[+-]?Integer>")
        parser.add_option("-p", type="float", dest="p_value", help="p-value threshold for chi-square test and t-test [default: %default]", default=1E-4, metavar="<Float>")
        parser.add_option("-c", "--coverage", type="int", dest="coverage", help="minimum coverage at mC position in DMR [default: %default]", default=5, metavar="<Integer>")
        parser.add_option("--filter", type="string", dest="filter", help="mC context to be called DMR. [default: %default]", default='CG', metavar="CG[, CHG[,CHH]")
        parser.add_option("-l", "--label", type="string", dest="label", help="camma-separated label strings for sample name in DMR output file [default: Sample1,Sample2]", metavar="String[,String]")
        parser.add_option("--neighbor-dist", type="int", dest="neighbor_dist", help="max distance between neighboring mC positions in DMR [default: %default]", default=100, metavar="<Integer>")
        parser.add_option("-m", "--multi-thread", type="int", dest="thread", help="number of threads [default: %default]", default=1, metavar="<Integer>")
        parser.add_option("--c-samples", type="int", dest="c_samples", help="minimum number of mCs for DMR [default: %default]", default=5, metavar="<Integer>")
        parser.add_option("--min-width", type="int", dest="min_width", help="minimum width of DMR [default: %default]", default=100, metavar="<Integer>")
        
        (opts, args) = parser.parse_args(argv)
        
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + str(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2
    
    if DEBUG:
        print ""
        print "smooth=",opts.smooth
        print "smooth_window=",opts.smooth_window
        print "target=",opts.target
        print "gtf=",opts.gtf
        print "assoc_dist=",opts.assoc_dist
        print "p_value=",opts.p_value
        print "coverage=",opts.coverage
        print "filter=",opts.filter
        print "label=",opts.label
        print "neighbor_dist=",opts.neighbor_dist
        print "thread=",opts.thread
        print "c_samples=",opts.c_samples
        print "min_width=",opts.min_width            
    
    ''' Validate Arguments '''
    if len(args) < 2:
        parser.print_usage()
        sys.exit(1)
    
    if opts.target != None and not os.path.isfile(opts.target):
        sys.stderr.write("Target region file %s (indicated by -t, --target) doesn't exist" % opts.target)
        sys.exit(1)
    
    ''' Prepare '''
    print os.getcwd()
    #tmpdir='./tmp_'+datetime.now().strftime('%Y%m%d%H%M%S')
    #os.mkdir(tmpdir)
    
    t_pool = ThreadPool(opts.thread)
    
    ''' Read Sample Files (Concurrent by sample)'''
    for m, sample_arg in enumerate(args):
        replist = sample_arg.split(',')
        for n, replicate_file in enumerate(replist):
            if not os.path.exists(replicate_file):
                sys.stderr.write('%d th replicate file of %d th sample (%s) doesn\'t exist' % (n+1,m+1,replicate_file))
                sys.exit(1)
        sample_list.append(SampleData(replist))
        
    for m, sample_data in enumerate(sample_list):
        t_pool.queueTask(preprocess_samples,sample_data,None)
    
    t_pool.joinAll()
    
    if DEBUG:
        print "print chromosome order"
        for chrom_name in sample_data.chrom_order:
            print chrom_name
    
    ''' debug purpose '''
    if DEBUG:
        for sample_data in sample_list:
            sample_data.output_debug_info()
    
    ''' Smoothing'''
    if opts.smooth:
        for sample_data in sample_list:
            smooth_data((sample_data,opts.smooth_window))
    
    ''' Filtering '''
    context_list = [field.strip() for field in opts.filter.split(",")]
    for sample_data in sample_list:
        sample_data.filter_by_context(context_list)
    
    ''' debug purpose '''
    if DEBUG:
        for sample_data in sample_list:
            sample_data.output_debug_info()
    
    ''' Detect DMRs (Concurrent by chromosome) '''
    gene_container = None
    if opts.gtf:
        gene_container = GeneParser.parse(opts.gtf)
        
    t_pool = ThreadPool(opts.thread)
    for m, sample_m in enumerate(sample_list):
        for n, sample_n in enumerate(sample_list[(m+1):]):
            dmr_dict = None
            if opts.target:
                if sample_m.get_number_of_replicates() < 3 and sample_n.get_number_of_replicates() < 3:
                    dmr_dict = SignificanceTest.test( t_pool, SignificanceTest.METHOD_TARGET_CHISQUARE, sample_m, sample_n, opts)
                else:
                    dmr_dict = SignificanceTest.test( t_pool, SignificanceTest.METHOD_TARGET_TTEST, sample_m, sample_n, opts)
            else:
                if sample_m.get_number_of_replicates() < 3 and sample_n.get_number_of_replicates() < 3:
                    dmr_dict = SignificanceTest.test( t_pool, SignificanceTest.METHOD_NOTARGET_CHISQUARE, sample_m, sample_n, opts)
                else:
                    dmr_dict = SignificanceTest.test( t_pool, SignificanceTest.METHOD_NOTARGET_TTEST, sample_m, sample_n, opts)
            
            assign_dmr_id(dmr_dict)
            if gene_container:
                associate_gene(gene_container, dmr_dict, 'diff-gene-%d-%d.txt'%(m+1,n+2), opts.assoc_dist, saved_argv)
            output_dmr(sample_m.chrom_order, dmr_dict, 'diff-dmr-%d-%d.txt'%(m+1,n+2), saved_argv)
    
    ''' remove temp directory'''
    #os.rmdir(tmpdir)

    print "finish bsf_diff"

if __name__ == "__main__":
    if DEBUG:
        pass
    if TESTRUN:
        import doctest
        doctest.testmod()
    sys.exit(main())
