import unittest
import os
from bsf_diff.gene import Gene, GeneContainer, GeneParser

class Test(unittest.TestCase):

    def setUp(self):
        self.testgtf = os.path.join(os.path.dirname(__file__),"test_data","test.gtf")
        self.errgtf = os.path.join(os.path.dirname(__file__),"test_data","error.gtf")
        print self.testgtf
        
    def tearDown(self):
        pass

    def testName(self):
        pass
    
    def testParse(self):
        geneContainer = GeneParser.parse(self.testgtf)
        geneChr1List = geneContainer.get_of_chrom('chr1')[0]
        print geneChr1List
        geneChr2List = geneContainer.get_of_chrom('chr2')[0]
        print geneChr2List
        self.assertEqual(geneChr1List.start, 66999825, 'start of NM_032291')
        self.assertEqual(geneChr1List.end, 67136702, 'end of NM_032291')
        self.assertEqual(geneChr2List.start, 201170985, 'start of NM_001100422')
        self.assertEqual(geneChr2List.end, 201334739, 'end of NM_001100422')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()