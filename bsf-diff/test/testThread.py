import unittest
import os
import sys
import time
from bsf_diff.threadpool import ThreadPool, ThreadPoolThread

def threadMethod(data):
    print data
    
class Test(unittest.TestCase):

    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def testName(self):
        pass
    
    def testParse(self):
        pool = ThreadPool(2)
        pool.queueTask(threadMethod, 'hogepiyo', None)
        pool.queueTask(threadMethod, 'test', None)
        pool.joinAll()



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()