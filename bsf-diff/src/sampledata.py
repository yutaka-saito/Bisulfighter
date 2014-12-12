'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''

from bsfdiffexception import BsfDiffException
from methylc import MethylCManager

class SampleData:

    def __init__(self, rep_file_list):
        self._num_rep=len(rep_file_list)
        self._rep_files=rep_file_list
        self._chrom_dict={} #key is chromosome name, value is MethylCManager
        self.chrom_order=[] #order of chromosome appearance in input file. used later
    
    def _get_chrom_methyl_c_manager(self, chrom):
        if chrom not in self._chrom_dict:
            self._chrom_dict[chrom] = MethylCManager(chrom,self._num_rep)
        return self._chrom_dict[chrom]
    
    def read_from_file(self):
        is_read_chrom_order = False # chromosome order has been read already or not
        for m, rep_file in enumerate(self._rep_files):
            line_num = 0
            current_chrom_name = ""
            for line in open(rep_file):
                line_num = line_num + 1
                line=line.strip()
                if line[0] == '#':
                    continue
                fields = line.split('\t')
                if len(fields) != 6:
                    raise BsfDiffException("The following line in %s is in wrong format,line number is %d"%(rep_file, line_num))
                
                ''' parse columns to variables '''
                chrom = fields[0]
                if not is_read_chrom_order and chrom != current_chrom_name:
                    self.chrom_order.append(chrom)
                    current_chrom_name = chrom
                pos = int(fields[1])
                is_plus_strand = True
                if fields[2] != '+':
                    is_plus_strand = False
                context = fields[3]
                methyl_rate = float(fields[4])
                coverage = fields[5]
                
                ''' put value to manager'''
                self._get_chrom_methyl_c_manager(chrom).add_replicate( m, pos, is_plus_strand, context, methyl_rate, coverage)
                
            is_read_chrom_order = True
    
    def filter_by_context(self, context_list):
        for manager in self._chrom_dict.values():
            manager.filter_by_context(context_list)
    
    def get_manager_of_chrom(self, chrom):
        if chrom in self._chrom_dict:
            return self._chrom_dict[chrom];
    
    def get_chromosome_list(self):
        return self._chrom_dict.keys()
    
    def get_number_of_replicates(self):
        return self._num_rep
    
    def output_debug_info(self):
        print 'output new SampleData'
        for chrom,manager in self._chrom_dict.items():
            print 'chrom='+chrom
            for methylc in manager:
                print methylc
    
    def __str__(self):
        return "self._rep_files %s"%self._rep_files