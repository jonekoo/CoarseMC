#-*- coding: utf-8 -*-
"""
Reorders PTGBCYL output by temperature. 

Usage example:
python reorder_output.py output-*.json

The script creates corresponding files sorted-0.json, sorted-1.json,
etc. Together they contain the same information as the input files of
the script but now a single file contains only the output at a single
temperature. The files are indexed in order of ascending temperature.

Note 1) ujson is used if it is available. Encoding is ~10 times faster
        with ujson as compared to the json module in Python 2.7.9.
Note 2) memory required is roughly the same as the size of the input
        files combined.
"""
import sys
try:
    # Encoding with ujson is around 10 times faster 
    import ujson as json
except ImportError:
    # than with this Python 2.7.9 module.
    import json

def main():
    fns = sys.argv[1:]
    all_results = {}
    for fn in fns:
        print "Processing file " + fn
        try:
            f = open(fn, 'r')
        except IOError:
            print __doc__
            return 1
        json_str = f.read()
        f.close()
        splitted = json_str.split("}\n{")
        try: 
            d = json.loads(splitted[0] + "}")
            if d["temperature"] in all_results:
                all_results[d["temperature"]].append(d)
            else:
                all_results[d["temperature"]] = [d]
        except ValueError:
            pass
        for item in splitted[1:-1]:
           try:
               d = json.loads("{" + item + "}")
           except ValueError:
               continue
           if d["temperature"] in all_results:
               all_results[d["temperature"]].append(d)
           else:
               all_results[d["temperature"]] = [d]
        try:
            d = json.loads("{" + splitted[-1])
            if d["temperature"] in all_results:
                all_results[d["temperature"]].append(d)
            else:
                all_results[d["temperature"]] = [d]
        except ValueError:
            pass

    for i, key in enumerate(sorted(all_results)):
        all_results[key].sort(key=lambda result: result["i_sweep"])
        f = open("sorted-" + str(i) + ".json", 'w')
        json.dump(all_results[key], f)
        f.close()



if __name__=='__main__':
    main()
