import os
import unittest
from bsf_diff.target_region import TargetRegion, TargetRegionContainer, TargetRegionParser

class Test(unittest.TestCase):
    
    def setUp(self):
        self.testbed = os.path.join(os.path.dirname(__file__),"test_data","test_normal.bed")
        
    def tearDown(self):
        pass
    
    def testName(self):
        trcontainer = TargetRegionParser.parse( self.testbed )
        target_region_list_chr1 = trcontainer.get_of_chrom("chr1")
        target_region_list_chr2 = trcontainer.get_of_chrom("chr2")
        
        self.assertEqual(target_region_list_chr1[0].chrom, 'chr1', 'chr of target1')
        self.assertEqual(target_region_list_chr1[0].start, 100, 'start of target1')
        self.assertEqual(target_region_list_chr1[0].end, 200, 'end of target1')
        self.assertEqual(target_region_list_chr1[0].name, 'target1', 'name of target1')
        
        self.assertEqual(target_region_list_chr1[1].chrom, 'chr1', 'chr of target2')
        self.assertEqual(target_region_list_chr1[1].start, 111, 'start of target2')
        self.assertEqual(target_region_list_chr1[1].end, 10000, 'end of target2')
        self.assertEqual(target_region_list_chr1[1].name, 'target2', 'name of target2')
        
        self.assertEqual(target_region_list_chr2[0].chrom, 'chr2', 'chr of target3')
        self.assertEqual(target_region_list_chr2[0].start, 100000, 'start of target3')
        self.assertEqual(target_region_list_chr2[0].end, 200000, 'end of target3')
        self.assertEqual(target_region_list_chr2[0].name, 'target3', 'name of target3')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()