#-*- coding: utf-8 -*-
import sys
import json

def main():
    filename = sys.argv[1]
    if len(sys.argv) > 2:
        try:
            outfile = open(sys.argv[2], 'w')
        except:
            outfile = None
    f = open(filename, 'r')
    i_configuration = 0
    rodgroup = {}
    rodgroup['type'] = 'rod'
    rodgroup['name'] = 'gb'
    rodgroup['description'] = ['x', 'y', 'z', 'ux', 'uy', 'uz']

    pointgroup = {}
    pointgroup['type'] = 'point'
    pointgroup['name'] = 'lj'
    pointgroup['description'] = ['x', 'y', 'z']
    for line in f:
        if line.startswith('#configuration'):
            i_configuration += 1
            rodgroup['coordinates'] = []
            pointgroup['coordinates'] = []
        elif line.startswith('#end_configuration'):
            particle_groups = [rodgroup]
            rodgroup['n'] = len(rodgroup['coordinates'])
            if (len(pointgroup['coordinates']) > 0):
                pointgroup['n'] = len(pointgroup['coordinates'])
                particle_groups.append(pointgroup)
            if outfile is not None:
                json.dump({'particle_groups': particle_groups}, outfile, sort_keys=False, indent=4, separators=(',', ': '))
            else:
                print json.dumps({'particle_groups': particle_groups}, sort_keys=False, indent=4, separators=(',', ': '))
        elif line.startswith('gb'):
            rodgroup['coordinates'].append(line.split()[1:])
        elif line.startswith('lj'):
            pointgroup['coordinates'].append(line.split()[1:4])
            
if __name__=='__main__':
    main()
