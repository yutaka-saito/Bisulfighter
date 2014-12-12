'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''
import numpy
from bsfdiffexception import BsfDiffException

class DMR:
    def __init__( self, chrom, start_pos, end_pos, strand, mc_rate_1, mc_rate_2, mc_stddev_1, mc_stddev_2, fold_change, p_value ):
        self.chrom = chrom
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.strand = strand
        self.mc_rate_1 = mc_rate_1
        self.mc_rate_2 = mc_rate_2
        self.mc_stddev_1 = mc_stddev_1
        self.mc_stddev_2 = mc_stddev_2
        self.fold_change = fold_change
        self.p_value = p_value
        self.dmr_id = None

class DMRCandidate:
    def __init__(self):
        self.methylc_m_list = list()
        self.methylc_n_list = list()
        self.p_value_list = list()
    
    def add_methylc(self, methylc_m, methylc_n, p_value):
        self.methylc_m_list.append( methylc_m )
        self.methylc_n_list.append( methylc_n )
        self.p_value_list.append( p_value )
        if methylc_m.get_average_methylrate > methylc_n.get_average_methylrate:
            is_first_higher = True
        else:
            is_first_higher = False
        if not hasattr(self, "_is_first_higher"):
            self._is_first_higher = is_first_higher
        elif self._is_first_higher != is_first_higher:
            raise BsfDiffException("new methyl_c's methyl rate is contradicted")
    
    def is_first_sample_is_higher(self):
        ''' True if first sample(m) has higher methylation rate, False otherwise'''
        return self._is_first_higher
    
    def get_number_of_site(self):
        return len(self.methylc_m_list)

    def get_start_and_end(self):
        start = None
        end = None
        for methyl_c in self.methylc_m_list:
            if not start or methyl_c.position < start:
                start = methyl_c.position
            if not end or end < methyl_c.position:
                end = methyl_c.position
        return (start,end)
    
    def get_length(self):
        (start, end) = self.get_start_and_end()
        return end - start + 1
    
    def get_strand(self):
        strand = None
        for methyl_c in self.methylc_m_list:
            if strand == None:
                strand = methyl_c.is_plus_strand()
            elif strand != methyl_c.is_plus_strand():
                return '.'
        if strand:
            return '+'
        else:
            return '-'

    def get_average_and_stddev(self):
        methylrate_m_list = list()
        
        for methylc_m in self.methylc_m_list:
            methylrate_m_list.append(methylc_m.get_average_methylrate())
            
        array_m = numpy.array(methylrate_m_list)
        
        methylrate_n_list = list()
        
        for methylc_n in self.methylc_n_list:
            methylrate_n_list.append(methylc_n.get_average_methylrate())
            
        array_n = numpy.array(methylrate_n_list)
        
        return ( numpy.average( array_m ), numpy.average( array_n ), numpy.std( array_m ), numpy.std( array_n ) )
    
    def get_average_pvalue(self):
        sum_p_value = 0
        for p_value in self.p_value_list:
            sum_p_value = sum_p_value + p_value
        
        return sum_p_value / len(self.p_value_list)  
        
class MethylC:
    ''' enum for methylation context '''
    CONTEXT_CG, CONTEXT_CHG, CONTEXT_CHH, CONTEXT_UNKNOWN = range(4)
    
    def __init__(self, num_replicate, position):
        self._rep_array = [None for i in xrange(num_replicate)]
        self.position = position
    
    @classmethod
    def context_str_to_enum(cls, str_context):
        if str_context == 'CG':
            return cls.CONTEXT_CG
        elif str_context == 'CHG':
            return cls.CONTEXT_CHG
        elif str_context == 'CHH':
            return cls.CONTEXT_CHH
        else:
            return cls.CONTEXT_UNKNOWN
    
    @classmethod
    def context_enum_to_str(cls, enum_context):
        if enum_context == cls.CONTEXT_CG:
            return 'CG'
        elif enum_context == cls.CONTEXT_CHG:
            return 'CHG'
        elif enum_context == cls.CONTEXT_CHH:
            return 'CHH'
        else:
            return 'Unknown'
    
    def add_replicate_value(self, replicate_index, is_plus_strand, context, methyl_rate, coverage):
        '''
        replicate_index: 0-based serial number of replicates
        is_plus_strand: True if this c is + strand, False otherwise
        context: methylation context. One of CG(MethylC.CONTEXT_CG), CHG(MethylC.CONTEXT_CHG), CHH(MethylC.CONTEXT_CHH)
        methyl_rate: rate of methylation at this postion
        coverage: read coverage at this position
        '''
        enumContext = MethylC.context_str_to_enum(context)
        self._rep_array[replicate_index] = [is_plus_strand, enumContext, float(methyl_rate), int(coverage)]
        
    def is_plus_strand(self):
        '''
        if plus strand, return True
        if minus strand, return False
        if unknown, return None
        '''
        if len(self._rep_array) == 0:
            return None
        return self._rep_array[0][0]
    
    def get_average_coverage(self):
        if len(self._rep_array) == 0:
            return None
        coverage_list = list()
        for rep in self._rep_array:
            coverage_list.append(rep[3])
        return numpy.average( numpy.array(coverage_list) )
    
    def get_average_methylrate(self):
        if len(self._rep_array) == 0:
            return None
        methylrate_list = list()
        for rep in self._rep_array:
            methylrate_list.append(rep[2])
        return numpy.average( numpy.array(methylrate_list) )
    
    def get_context(self):
        if len(self._rep_array) == 0:
            return None
        return MethylC.context_enum_to_str(self._rep_array[0][1])
    
    def get_methylc_nonmethylc(self):
        '''
        return (number of methyl-C, number of non-methyl-C)
        '''
        sum_met_c = 0
        sum_nonmet_c = 0
        for rep in self._rep_array:
            methyl_rate = rep[2]
            coverage = rep[3]
            num_met_c = int(round(coverage * methyl_rate))
            num_nonmet_c = coverage - num_met_c
            
            sum_met_c += num_met_c
            sum_nonmet_c += num_nonmet_c
            
        return (sum_met_c, sum_nonmet_c)
    
    def get_methylrate_list(self):
        methylrate_list = list()
        
        for rep in self._rep_array:
            if rep[2]:
                methylrate_list.append(rep[2])
        
        return methylrate_list
    
    def get_methylrate(self,rep_num):
        if rep_num < len(self._rep_array):
            rep = self._rep_array[rep_num]
            return rep[2]
        else:
            return None
    
    def get_coverage_list(self):
        coverage_list = list()
        
        for rep in self._rep_array:
            if rep[3]:
                coverage_list.append(rep[3])
        
        return coverage_list
    
    def get_coverage(self,rep_num):
        if rep_num < len(self._rep_array):
            rep = self._rep_array[rep_num]
            return rep[3]
        else:
            return None
    
    def update_methylrate(self,rep_num,methylrate):
        self._rep_array[rep_num][2]=methylrate
    
    def __str__(self):
        returnStr = str(' position=%d,'%self.position)
        
        for i, values in enumerate(self._rep_array):
            returnStr = returnStr+'(rep='+str(i)+', '
            if values:
                returnStr = returnStr+str('isplus=%s, context=%s, methyl_rate=%f, coverage=%d'%(values[0], MethylC.context_enum_to_str(values[1]), values[2], values[3]))+')'
            
            else:
                returnStr = returnStr+'None)'
        
        return str(returnStr)

class MethylCManager:
    ''' MethylCManager contains methylation state of A chromosome
    behave like iterator for all methylation state(=MethylC)
    '''
    BASES_PER_UNIT = 1000
    
    def __init__(self, chrom, num_replicate):
        self.chrom = chrom
        self.num_replicate = num_replicate
        self._methylc_list = []
        self._sorted = False
    
    ''' private methods '''
    def _add_methyl_c(self, position, methylc):
        list_index = self._get_array_index(position)
        self._methylc_list[list_index].append(methylc)
    
    def _get_array_index(self, position):
        list_index = position / MethylCManager.BASES_PER_UNIT
        for i in range(len(self._methylc_list),list_index+1):
            self._methylc_list.append([])
        return list_index
    
    ''' public methods '''
    def add_replicate(self, replicate_index, position, is_plus_strand, context, methyl_rate, coverage):
        self._sorted = False
        
        methyl_c = self.get_methyl_c(position)
        if not methyl_c:
            methyl_c = MethylC(self.num_replicate, position)
            self._add_methyl_c(position, methyl_c)
        methyl_c.add_replicate_value(replicate_index, is_plus_strand, context, methyl_rate, coverage)
    
    def filter_by_context(self, context_list):
        for eachlist in self._methylc_list:
            for i, each_methylc in enumerate(eachlist):
                to_be_filtered = True
                context = each_methylc.get_context()
                for filt_context in context_list:
                    if context == filt_context:
                        to_be_filtered = False
                        break
                if to_be_filtered:
                    eachlist[i] = 'to_be_deleted'
            for i in range(0,eachlist.count('to_be_deleted')):
                eachlist.remove('to_be_deleted')
    
    def get_methyl_c(self, position):
        array_index = self._get_array_index(position)
        
        ''' Out of range '''
        if len(self._methylc_list) < (array_index + 1):
            return None
        
        for methylc_cand in self._methylc_list[array_index]:
            if methylc_cand.position == position:
                return methylc_cand
        
        return None
    
    '''Iterator-related methods'''
    def __iter__(self):
        if not self._sorted:
            for i, methylcs in enumerate(self._methylc_list):
                self._methylc_list[i] = sorted(methylcs, key=lambda methylc: methylc.position)
            self._sorted = True
        
        self._arrayIndexForIter=0
        self._IndexInArrayIndexForIter=0
        return self;
    
    def next(self):
        isFound = False
        returnObj = None
        while not isFound:
            if len(self._methylc_list) <= self._arrayIndexForIter:
                raise StopIteration
            
            list = self._methylc_list[self._arrayIndexForIter]
            
            if len(list) <= self._IndexInArrayIndexForIter:
                self._arrayIndexForIter = self._arrayIndexForIter + 1
                self._IndexInArrayIndexForIter = 0
                continue
            else:
                isFound = True
                returnObj = self._methylc_list[self._arrayIndexForIter][self._IndexInArrayIndexForIter]
                self._IndexInArrayIndexForIter = self._IndexInArrayIndexForIter + 1
        
        return returnObj
            
