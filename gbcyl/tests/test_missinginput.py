#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Test that simulation is aborted correctly when particle_groups is missing
from input.
"""
import json
import subprocess
import sys
import os.path


def main():
    # Load valid input JSON.
    with open('test_missinginput.json', 'r') as f:
        d = json.load(f)
        # Remove particle_groups from JSON.
        del d['particle_groups']
    # Save JSON.
    with open('test_missinginput-0.json', 'w') as f:
        json.dump(d, f)
    # Run program with input. Catch output.
    args = ["mpirun", "-np", "1", os.path.join("..", "src", "ptgbcyl"), "-i",
            "test_missinginput-_I_.json"]
    proc = subprocess.Popen(args, stderr=subprocess.PIPE)
    for line in proc.stderr:
        if "error" in line.lower() and "particle_groups" in line:
            sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
