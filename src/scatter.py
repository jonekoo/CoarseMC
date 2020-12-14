#-*- coding: utf-8 -*-
"""
Scatters the variables prefixed with _SCATTER_ in the master file
across a set of files.
"""
import sys
import json



def main():
    prefix = '_SCATTER_'
    with open(sys.argv[1], 'r') as f:
        master = json.load(f)
    fn = sys.argv[2]
    n = int(sys.argv[3])
    idstr = '_I_'
    # load data to scatter
    # walk through keys
    for id in range(n):
        with open(fn.replace(idstr, str(id)), 'r') as f:
            d = json.load(f)
            for key in master:
                if key.startswith(prefix):
                # pick the right value for the current replica and
                    # substitute
                    d[key.replace(prefix, '')] = master[key][id]
        with open(fn.replace(idstr, str(id)), 'w') as f:
            json.dump(d, f)


if __name__ == '__main__':
    main()
