#-*- coding: utf-8 -*-
import sys
import json

descriptions = []

def main():
    f = open(sys.argv[1], 'r')
    d = json.load(f, object_hook=delete_others)
    walk(d, 0, '') #, descriptions)


def delete_others(d):
    for key, value in d.items():
        if isinstance(value, dict) or isinstance(value, list):
            pass
        elif " description" not in key:
            del d[key]
    return d
                
def parse_pairs(pairs):
    for p in pairs:
        if " description" in p[0]:
            descriptions.append((p[0].split()[0], p[1]))
    

def walk(d, indent, parameter):
    if isinstance(d, dict):
        print "\\begin{itemize}"
        indent = indent + 4
        for key, value in d.items():
            parameter = key.split()[0]
            print " " * indent + "\\item[" + parameter.replace('_', '\\_') + ": ]",
            walk(value, indent, parameter)
        print "\\end{itemize}"
    elif isinstance(d, list):
        for item in d:
            walk(item, indent, parameter)
    elif isinstance(d, unicode):    
        print d.replace(parameter, parameter.replace('_', '\\_'))
                

if __name__=='__main__':
    main()

