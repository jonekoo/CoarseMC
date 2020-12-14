#-*- coding: iso-8859-1 -*-
import numpy
import replaceparameter
import sys
import subprocess

# read temperatures
betas=numpy.loadtxt('beta.dat')
temperatures=[1.0/beta for beta in betas]
temperatures.sort()
ntemperatures=len(temperatures)
ntasks=2
for i in range(ntemperatures-ntasks, -1, -1):
    newtemps=temperatures[i:i+ntasks]
    for j in range(ntasks):
        p=subprocess.call(["cp", "restartparameters."+str(j), "inputparameters."+str(j)])
        p=subprocess.call(["cp", "restartconfiguration."+str(j), "inputconfiguration."+str(j)])
        sys.argv[1:]=['inputparameters.'+str(j), 'temperature', str(newtemps[j])]
        replaceparameter.main() 
        sys.argv[1:]=['inputparameters.'+str(j), 'i_sweep', str(0)]
        replaceparameter.main() 
        sys.argv[1:]=['inputparameters.'+str(j), 'n_production_sweeps', str(0)]
        replaceparameter.main() 
        sys.argv[1:]=['inputparameters.'+str(j), 'n_equilibration_sweeps', str(1)]
        replaceparameter.main() 

    p=subprocess.call(["mpirun", "-np", str(ntasks), "../src/ptgbcyl"])

    if i>0:
        p=subprocess.call(["cp", "restartparameters."+str(ntasks-1), "restartparameters."+str(i+ntasks-1)])
        p=subprocess.call(["cp", "restartconfiguration."+str(ntasks-1), "restartconfiguration."+str(i+ntasks-1)])
