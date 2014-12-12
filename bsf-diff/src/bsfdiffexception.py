'''
bsf_diff: Defferentially Methylated Region Detection Tool For Bisulphite Data
Copyright 2013 Computational Biology Research Center, AIST. All rights reserved.
'''

class BsfDiffException(Exception):

    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return str(self.value)
        