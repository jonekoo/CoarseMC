#-*- coding: utf-8 -*-
"""
Usage:
python update_json.py master.json slave1.json [slave2.json ...]

Updates the JSON in files slave1.json slave2.json etc. with the JSON in
master.json.
"""
import sys
import json

def main():
    # Load "master" json used for updating
    fn_master = sys.argv[1]
    f = open(fn_master, 'r')
    d_master = json.load(f)
    f.close()
    # Load "slave" json to be updated
    fn_slaves = sys.argv[2:]
    for fn_slave in fn_slaves:
        f = open(fn_slave, 'r')
        d_slave = json.load(f)
        f.close()

        d_slave.update(d_master)

        # Overwrite "slave" file.
        f = open(fn_slave, 'w')
        json.dump(d_slave, f)
        f.close()

if __name__=='__main__':
    main()
