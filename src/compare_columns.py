#-*- coding: utf-8 -*-
"""
Prints differences of columns in two files. If only one column number is
given, it is used for both files.

Usage
python compare_columns.py filename1 filename2 #col [#col2]


"""
import sys
import numpy



def main():
    try:
        fn1 = sys.argv[1]
        fn2 = sys.argv[2]
        col1 = int(sys.argv[3])
        if len(sys.argv) > 4:
            col2 = int(sys.argv[4])
        else:
            col2 = col1
        data1 = numpy.loadtxt(fn1, usecols=[col1])
        data2 = numpy.loadtxt(fn2, usecols=[col2])
    except IOError:
        print __doc__
        return 1
    diff = data1 - data2
    for d1, d2, delta in zip(data1, data2, diff):
        print d1, d2, delta, delta / d1
    
    

if __name__=='__main__':
    main()
