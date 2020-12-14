# -*- coding: utf-8 -*-
import sys
import random

f = open(sys.argv[1], 'r')
n_gb = int(sys.argv[2])
n_lj = int(sys.argv[3])


lj_indices = random.sample(xrange(n_gb), n_lj)

lines = f.readlines()
i = 0
# three first lines are not molecules. 
for i in lj_indices:
    lines[i+3] = lines[i+3].replace('gb', 'lj')

for line in lines:
    print line.rstrip()

