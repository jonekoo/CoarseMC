#-*- coding: utf-8 -*-
import sys

def main():
    filename = sys.argv[1]
    f = open(filename, 'r')
    i_configuration = 0
    rodgroup = {}
    rodgroup['type'] = 'rod'
    rodgroup['name'] = 'gb'
    rodgroup['description'] = ['x', 'y', 'z', 'ux', 'uy', 'uz']

    pointgroup = {}
    pointgroup['type'] = 'point'
    pointgroup['name'] = 'gb'
    pointgroup['description'] = ['x', 'y', 'z', 'ux', 'uy', 'uz']
    while True:
        line = f.readline()
        if line.startswith('#configuration'):
            i_configuration += 1
            rodgroup['coordinates'] = []
            pointgroup['coordinates'] = []
        elif line.startswith('#end_configuration'):
            particle_groups = [rodgroup]
            if (len(pointgroup['coordinates']) > 0):
                particle_groups.append(pointgroup)
            json.dumps({'particle_groups': particle_groups})
        elif line.startswith('gb'):
            rodgroup['coordinates'].append(line.split()[1:])
        elif line.startswith('lj'):
            pointgroup['coordinates'].append(line.split()[1:4])

if __name__=='__main__':
    main()
