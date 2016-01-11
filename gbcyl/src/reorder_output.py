#-*- coding: utf-8 -*-
"""
Reorders PTGBCYL output by temperature. 

Usage example:
python reorder_output.py output-*.json

The script creates corresponding files sorted-0.json, sorted-1.json,
etc. Together they contain the same information as the input files of
the script but now a single file contains only the output at a single
temperature. The files are indexed in order of ascending temperature.

Note: memory required is about the same as the input files combined.
"""
import sys
import json

def main():
    fns = sys.argv[1:]
    all_results = {}
    for fn in fns:
        print "Processing file " + fn
        f = open(fn, 'r')
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
