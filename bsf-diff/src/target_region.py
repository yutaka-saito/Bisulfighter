'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''

import sys

class TargetRegion:
    
    def __init__(self, chrom, start, end, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
    
    def __str__(self):
        return "chrom=%s,start=%d,end=%d,name=%s" % (self.chrom, self.start, self.end, self.name)

class TargetRegionContainer:
    
    def __init__(self):
        self.target_region_dict = {} # key is chrom, value is list of target regions
        pass
    
    def add(self, target_region):
        if target_region.chrom in self.target_region_dict:
            self.target_region_dict[target_region.chrom].append(target_region)
        else:
            self.target_region_dict[target_region.chrom] = [target_region]
    
    def get_of_chrom(self, chrom):
        if chrom in self.target_region_dict:
            return self.target_region_dict[chrom]
        else:
            return list()

class TargetRegionParser:
    
    @classmethod
    def parse( self, target_file_name ):
        target_region_container = TargetRegionContainer()
        f = open(target_file_name,'r')
        try:
            line_num = 0
            for line in f:
                line_num += 1
                s_chrom, s_start, s_end, s_name = (line.strip().split('\t')+[None]*4)[:4]
                
                if s_name == None:
                    sys.stderr.write('line %d of BED file is wrong(name is absent)' % line_num)
                    f.close()
                    sys.exit(1)
                
                try:
                    target_region_container.add( TargetRegion( s_chrom, int(s_start), int(s_end), s_name ) )
                
                except ValueError:
                    sys.stderr.write('line %d of BED file is wrong' % line_num)
                    f.close()
                    sys.exit(1)
        
        finally:
            f.close()
        
        return target_region_container
        
        