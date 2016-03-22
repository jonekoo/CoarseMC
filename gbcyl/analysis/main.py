#-*- coding: utf-8 -*-
"""
This program computes the orientation parameter and director for a group
of particles. The program feeds coordinate data from JSON to the Fortran
routines in the library wrapper.so, which then performs the actual
computation.

Usage example: python main.py output-0.json

Note: you may need to define PYTHONPATH to point to the location of 
wrapper.so
"""
import wrapper
import json
import numpy
import sys
import matplotlib.pyplot as plt

def main():
    parameters = []
    try:
        f = open(sys.argv[1], 'r')
    except IOError:
        print __doc__
    else:
        with f:
            for line in f:
                d = json.loads(line)
                n = d["particle_groups"][0]["n"]
                coordinates = numpy.array(
                    d["particle_groups"][0]["coordinates"],
                    order='F').transpose()
                p2, director = wrapper.wrap_particledat(coordinates, n)
                print p2, director[0], director[1], director[2]
                #parameters.append(p2, director[0], director[1], director[2])

        
            
if __name__ == '__main__':
    main()
