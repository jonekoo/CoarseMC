#-*- coding: utf-8 -*-
"""
Replaces randomly selected particles from the first particle_group with
point-particles. A new group called 'lj' is created for the particles.

Usage: python replace_with_points.py n file1.json [file2.json ...]
"""
import sys
import json
import random

def main():
    k = int(sys.argv[1])
    pointgroup = {"name": "lj", "type": "point",
                  "description": ["x", "y", "z"], "n": k}
    for fn in sys.argv[2:]:
        with open(fn, 'r') as f:
            d = json.load(f)
        js = random.sample(range(len(d["particle_groups"][0]["coordinates"])),
                           k)
        coordinates = []
        for j in js:
            coordinates.append(d["particle_groups"][0]["coordinates"].pop(j))
        d["particle_groups"][0]["n"] = d["particle_groups"][0]["n"] - k
        pointgroup["coordinates"] = [item[:3] for item in coordinates]
        d["particle_groups"].append(pointgroup)
        with open(fn, 'w') as f:
            json.dump(d, f)

if __name__ == '__main__':
    main()
