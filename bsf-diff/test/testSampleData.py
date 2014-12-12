import unittest
import os
from bsf_diff.sampledata import SampleData
from bsf_diff.bsfdiffexception import BsfDiffException

def normal():
    pass

class Test(unittest.TestCase):

    def setUp(self):
        self.tsvwrong_numcolumn = os.path.join(os.path.dirname(__file__),"test_data","error_call_1.tsv")
            
    def tearDown(self):
        pass
    
    def testError1(self):
        sampleData = SampleData([self.tsvwrong_numcolumn])

        self.assertRaises(BsfDiffException, sampleData.read_from_file)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()