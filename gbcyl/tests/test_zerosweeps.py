#!/usr/bin/env python 
#-*- coding: utf-8 -*-
import subprocess
import json
import os.path
import sys
import numpy

def main():
    program = os.path.realpath("ptgbcyl")
    #inputdir = os.path.join("tests",
    #                        "zero_sweeps")
    inputdir = ""
    ntasks = 1
    energies = collect_energies(ntasks, inputdir)
    compare_energies(ntasks, inputdir, energies, tolerance = 1e-8)

def cleanup(ntasks, inputdir):
    for i in range(ntasks):
        subprocess.call(["rm", "-rf", str(i)], cwd=inputdir)
        

def initialize(ntasks, inputdir):
    # Create temporary directories for the input files.
    for i in range(ntasks):
        subprocess.call(["mkdir", "-p", str(i)], cwd=inputdir)
        
    # Copy input files to the directories.
    for i in range(ntasks):
        subprocess.call(["cp", os.path.join(inputdir,
                        "inputconfiguration-inputparameters-" + str(i) + \
                                            ".json"), os.path.join(str(i),
                                            "input-0.json")], cwd=inputdir)


def run_program(ntasks, program, inputdir):
    # Run program for each input files
    for i in range(ntasks):
        subprocess.call(["mpirun", "-np", "1", program],
                        cwd=os.path.join(inputdir, str(i)))
                            
def collect_energies(ntasks, inputdir):
    energies = []
    for i in range(ntasks):
        f = open("zerosweepsrestart-" + str(i) + ".json")
        d = json.load(f)
        energies.append(d["total_energy"])
        f.close()
    return energies

def compare_energies(ntasks, inputdir, energies, tolerance):
    reference = numpy.loadtxt(os.path.join(os.path.dirname(__file__),
                                           "total_energies-reference.txt"))
    # accept if differences are below tolerance
    accept = True
    for i, t in enumerate(zip(energies, reference[:ntasks])):
        e, r = t
        relative_error = abs((e - r) / r)
        if relative_error > tolerance:
             accept = False
             break
        else:
            print "Relative error = ", relative_error
    if accept:
        print "All energies are OK!"
    else:
        print "FAILED: Relative error in total_energy for task " + str(i) + \
            " = " + str(relative_error) + " > " + str(tolerance)
        sys.exit(1)


if __name__ == '__main__':
   main()
