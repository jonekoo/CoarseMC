#-*- coding: utf-8 -*-
import sys
import json

def main():
    # Load "master" json used for updating
    fn_master = sys.argv[1]
    f = open(fn_master, 'r')
    d_master = json.load(f)
    f.close()
    # Load "slave" json to be updated
    fn_slave = sys.argv[2]
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
