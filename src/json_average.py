#-*- coding: utf-8 -*-
"""
This script computes the average total energy for CoarseMC output sorted
by temperature. Before using this, one should first process the files with 
reorder_output.py.

Usage example:
python json_average.py sorted-0.json

"""
import sys
import json
#import blocked_error
import numpy

def main():
    filename = sys.argv[1]
    try:
        f = open(filename, 'r')
    except IOError:
        print __doc__
        return 1
    d = json.load(f)
    f.close()
    n = 0
    temperature = d[0]["temperature"]
    for g in d[0]["particle_groups"]:
        n += len(g["coordinates"])
    print filename, ": ", "T= ", d[0]["temperature"], "N= ", n, "<E_tot>/N= ", 
    data = []
    for item in d:
        data.append(item["total_energy"])
        if item["temperature"] != temperature:
            print 'ERROR: file contains several temperatures. See usage.'
            return
    print numpy.average(data) / n
    
        
if __name__=='__main__':
    main()
