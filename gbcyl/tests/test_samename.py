#!/usr/bin/env python 
#-*- coding: utf-8 -*-

"""
Check that the program fails correctly if two groups have the same name.
"""
import subprocess
import os.path
import sys

def main():
    runcmd = os.environ["RUN_CMD"].replace('"', '')
    args = runcmd.split() + ["1", os.path.join("..", "src", "ptgbcyl")]
    p = subprocess.Popen(args, stderr=subprocess.PIPE)
    for line in p.stderr:
        if 'error' in line.lower() and 'same name' in line.lower():
            sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
