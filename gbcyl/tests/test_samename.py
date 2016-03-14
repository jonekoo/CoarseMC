#-*- coding: utf-8 -*-
"""
Check that the program fails if two groups have the same name.
"""
import subprocess
import os.path
import sys

def main():
    if subprocess.call(["mpirun", "-np", "1",
                        os.path.join("..", "src", "ptgbcyl")]) == 1:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
