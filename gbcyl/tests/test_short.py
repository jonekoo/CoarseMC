#-*- coding: utf-8 -*-
"""
Tests running a short simulation 
"""
import subprocess
import json
import os.path
import sys
import numpy

def main():
    #inputdir = os.path.join("tests", "short")
    
    ntasks = 2
    tolerance = 1e-8
    for i in range(ntasks):
        reference_path = os.path.join(os.path.dirname(__file__), 
                                      "shortreference-" + str(i) + ".json")
        output_path = os.path.join("shortrestart-"+ str(i) + ".json")
        with open(reference_path) as f:
            reference = json.load(f)
        with open(output_path) as f:
            output = json.load(f)
        compare(output, reference, tolerance)

def compare(output, reference, tolerance):
    o = output["enthalpy"]
    r = reference["enthalpy"]
    if abs((o - r) / r) > tolerance:
       print "FAILED: Relative error in enthalpy = ", (o - r) / r  
       sys.exit(1)
       
if __name__ == '__main__':
   main()
