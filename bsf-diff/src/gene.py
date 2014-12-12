'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''
from bsfdiffexception import BsfDiffException

class Gene:
    
    def __init__(self, chrom, gene_id, start, end, strand):
        self.chrom = chrom
        self.gene_id = gene_id
        self.start = start
        self.end = end
        if "+" == strand:
            self.is_plus_strand = True
        else:
            self.is_plus_strand = False
    
    def __str__(self):
        return "chrom=%s,geneId=%s,start=%d,end=%d,strand=%s" % (self.chrom,self.gene_id,self.start,self.end,self.is_plus_strand)

class GeneContainer:

    def __init__(self):
        self.gene_dict = {} # key is chrom, value is list of genes
        pass
    
    def add(self, gene):
        if gene.chrom in self.gene_dict:
            self.gene_dict[gene.chrom].append(gene)
        else:
            self.gene_dict[gene.chrom] = [gene]
    
    def get_of_chrom(self, chrom):
        return self.gene_dict[chrom]
    
class GeneParser:

    @classmethod
    def parse( self, gtf_file_name ):
        gene_container = GeneContainer()
        f = open(gtf_file_name,'r')
        gene_dict = {} # key is chrom+gene_id, value is Gene instance
        try:
            line_num = 0
            for line in f:
                line_num = line_num + 1
                s_seqname, s_source, s_feature, s_start, s_end, s_score, s_strand, s_frame, s_attributes = (line.strip().split('\t')+[None]*9)[:9]
                
                i_start = int(s_start)
                i_end = int(s_end)
                
                for attr in s_attributes.split(";"):
                    if len(attr) == 0:
                        continue
                    attr_name, attr_value = attr.strip().split()
                    if attr_name == 'gene_id':
                        gene_id = attr_value.replace("\"","")
                
                if (s_seqname+gene_id) in gene_dict:
                    gene_in_dict = gene_dict[s_seqname+gene_id]
                    if s_seqname != gene_in_dict.chrom:
                        raise BsfDiffException("seqname of %s is inconsistent with previous line in line %d", gene_id, line_num)
                    if i_start < gene_in_dict.start:
                        gene_in_dict.start = i_start
                    if gene_in_dict.end < i_end:
                        gene_in_dict.end = i_end
                        
                else:
                    gene_dict[s_seqname+gene_id] = Gene(s_seqname, gene_id, i_start, i_end, s_strand)
        
        finally:
            f.close()
        
        for each_gene in gene_dict.itervalues():
            gene_container.add(each_gene)
        
        return gene_container
        
        