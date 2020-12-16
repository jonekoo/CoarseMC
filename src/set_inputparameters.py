#-*- coding: utf-8 -*-
import sys
import replaceparameter
import os.path

parameter = sys.argv[1]
value = sys.argv[2]

for filename in  os.listdir("."):
    if filename.split(".")[0] == "inputparameters" and \
    unicode(filename.split(".")[1]).isnumeric():
        replaceparameter.replace(filename, parameter, value)

