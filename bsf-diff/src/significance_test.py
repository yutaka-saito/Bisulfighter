'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''

from target_region import TargetRegionParser
from bsfdiffexception import BsfDiffException
from methylc import DMRCandidate, DMR
from rpy import r

class SignificanceTest:
    
    METHOD_NOTARGET_CHISQUARE, METHOD_NOTARGET_TTEST, METHOD_TARGET_CHISQUARE, METHOD_TARGET_TTEST = range(4)
    
    @classmethod
    def test(cls, threadPool, method, sample_m, sample_n, options ):
        dmr_dict = {}
        
        target_region_container = None
        if options.target != None:
            target_region_container = TargetRegionParser.parse( options.target )
        
        for chrom in sample_m.get_chromosome_list():
            manager_m = sample_m.get_manager_of_chrom(chrom)
            manager_n = sample_n.get_manager_of_chrom(chrom)
            if method == cls.METHOD_TARGET_CHISQUARE or method == cls.METHOD_TARGET_TTEST:
                threadPool.queueTask(cls.test_by_target,(chrom, method, target_region_container.get_of_chrom(chrom), dmr_dict, manager_m, manager_n, options.target, options),None)
            
            elif method == cls.METHOD_NOTARGET_CHISQUARE or method == cls.METHOD_NOTARGET_TTEST:
                threadPool.queueTask(cls.test_by_detect_region,(chrom, method, dmr_dict, manager_m, manager_n, options),None)
            
        threadPool.joinAll()
        
        return dmr_dict
    
    @classmethod
    def test_by_target(cls, (chrom, method, target_list, dmr_dict, manager_m, manager_n, target, options)):
        
        dmr_list = list()
        for target_region in target_list:
            
            is_dmr = True
            candidate_region = DMRCandidate()
            for methylc_m in manager_m:
                if methylc_m.position < target_region.start:
                    continue
                if target_region.end < methylc_m.position:
                    break
                
                methylc_n = manager_n.get_methyl_c( methylc_m.position )
                
                if methylc_n == None:
                    continue
                
                ''' Same sample must be higher than other '''
                if methylc_m.get_average_methylrate > methylc_n.get_average_methylrate:
                    is_first_higher = True
                else:
                    is_first_higher = False
                if candidate_region.get_number_of_site() != 0 and candidate_region.is_first_sample_is_higher() != is_first_higher:
                    is_dmr = False
                    break
                
                ''' Check coverage '''
                if methylc_m.get_average_coverage() < options.coverage or methylc_n.get_average_coverage() < options.coverage:
                    is_dmr = False
                    break
                
                ''' Check P-value '''
                p_value = None
                if method == cls.METHOD_TARGET_CHISQUARE:
                    p_value = cls._chisquare_for_csite(methylc_m, methylc_n);
                elif method == cls.METHOD_TARGET_TTEST:
                    p_value = cls._ttest_for_csite(methylc_m, methylc_n);
                
                if p_value <= options.p_value:
                    candidate_region.add_methylc(methylc_m, methylc_n, p_value)
            
            if is_dmr:
                ''' check the number of mC-site '''
                if candidate_region.get_number_of_site() < options.c_samples:
                    continue
                strand = candidate_region.get_strand()
                (rate_avg_m, rate_avg_n, rate_stddev_m, rate_stddev_n) = candidate_region.get_average_and_stddev()
                fold_change = rate_avg_m / rate_avg_n
                avg_p_value = candidate_region.get_average_pvalue()
                dmr_list.append(DMR(chrom,target_region.start,target_region.end,strand,rate_avg_m, rate_avg_n, rate_stddev_m, rate_stddev_n, fold_change, avg_p_value))
                
        dmr_dict[chrom] = dmr_list
    
    @classmethod
    def test_by_detect_region(cls, (chrom, method, dmr_dict, manager_m, manager_n, options)):
        
        dmr_list = list()
        current_candidate_region = None
        
        last_methylc_position = 0
        
        for methylc_m in manager_m:
            
            methylc_n = manager_n.get_methyl_c( methylc_m.position )
            
            if methylc_n == None:
                continue
            
            p_value = None;
            if method == cls.METHOD_NOTARGET_CHISQUARE:
                p_value = cls._chisquare_for_csite(methylc_m, methylc_n);
            elif method == cls.METHOD_NOTARGET_TTEST:
                p_value = cls._ttest_for_csite(methylc_m, methylc_n);
            else:
                raise BsfDiffException("No such test method %s", method)
            
            if current_candidate_region != None and options.neighbor_dist < (methylc_m.position - last_methylc_position - 1) :
                ''' if distance between this site and last site is larger than neighbor-dist, finish current candidate '''
                cls.validate_dmr(dmr_list, chrom, current_candidate_region, options)
                current_candidate_region = None
            
            last_methylc_position = methylc_m.position
            
            if current_candidate_region == None:
                ''' search for new region '''
                if p_value <= options.p_value:
                    current_candidate_region = DMRCandidate()
                    current_candidate_region.add_methylc( methylc_m, methylc_n, p_value )
            
            else:
                ''' try to recruit new mC to current candidate region '''
                if methylc_m.get_average_methylrate > methylc_n.get_average_methylrate:
                    is_first_higher = True
                else:
                    is_first_higher = False
                if p_value <= options.p_value and current_candidate_region.is_first_sample_is_higher() == is_first_higher:
                    current_candidate_region.add_methylc(methylc_m, methylc_n, p_value)
                else:
                    cls.validate_dmr(dmr_list, chrom, current_candidate_region, options)
                    current_candidate_region = None
        
        if current_candidate_region:
            cls.validate_dmr(dmr_list, chrom, current_candidate_region, options)
        
        dmr_dict[chrom] = dmr_list
    
    @classmethod
    def validate_dmr(cls, dmr_list, chrom, candidate_region, options):
        if options.c_samples <= candidate_region.get_number_of_site() \
           and options.min_width <= candidate_region.get_length():
            (start_pos, end_pos) = candidate_region.get_start_and_end()
            strand = candidate_region.get_strand()
            (rate_avg_m, rate_avg_n, rate_stddev_m, rate_stddev_n) = candidate_region.get_average_and_stddev()
            fold_change = None
            if rate_avg_n == 0:
                fold_change = "Inf"
            else:
                fold_change = rate_avg_m / rate_avg_n
            avg_p_value = candidate_region.get_average_pvalue()
            dmr_list.append(DMR(chrom,start_pos,end_pos,strand,rate_avg_m, rate_avg_n, rate_stddev_m, rate_stddev_n, fold_change, avg_p_value))
    
    @classmethod
    def _chisquare_for_csite(cls, methylc_m, methylc_n):
        #print "chisquare:m=",methylc_m.get_methylc_nonmethylc(),", n=",methylc_n.get_methylc_nonmethylc()
        matrix=r.cbind( methylc_m.get_methylc_nonmethylc(), methylc_n.get_methylc_nonmethylc() )
        return r.chisq_test(matrix,correct=False)['p.value']
    
    @classmethod
    def _ttest_for_csite(cls, methylc_m, methylc_n):
        
        return r.t_test(methylc_m.get_methylrate_list(),methylc_n.get_methylrate_list())
    