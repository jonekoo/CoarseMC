# -*- coding: utf-8 -*-
"""
Creates a geometric progression of temperatures and writes in JSON with
the key "pt_temperatures".
"""
import sys
import math
import json

def main():
    if len(sys.argv) < 4:
        print "usage: python geometric_progression.py lowest_temperature highest_temperature number_of_temperatures"
        sys.exit(1)
    else:
        # Give first and last number and the total number of points
        # including the ends.
        first = float(sys.argv[1])
        last = float(sys.argv[2])
        n = int(sys.argv[3])

        # Ratio of two consecutive numbers in the series
        ratio = math.pow(last/first, 1.0/(n-1))
        
        series = [ratio**i*first for i in range(n)]
        print json.dumps({'pt_temperatures': series})

if __name__ == '__main__':
    main()
