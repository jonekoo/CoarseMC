#-*- coding: utf-8 -*-
import wrapper
import json
import numpy
import sys
import matplotlib.pyplot as plt

def main():
    parameters = []
    with open(sys.argv[1], 'r') as f:
        for line in f:
            d = json.loads(line)
            n = d["particle_groups"][0]["n"]
            coordinates = numpy.array(d["particle_groups"][0]["coordinates"],
                                      order='F').transpose()
            p2, director = wrapper.wrap_particledat(coordinates, n)
            print p2, director
            parameters.append(p2)

            
if __name__ == '__main__':
    main()
