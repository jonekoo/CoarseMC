#!/bin/bash

# :NOTE: The $PATH environment variable must point to the directory where the 
# blocked_error exists or the blocked_error program must be in the same folder
# where this script is run from. 

rdf_file=radial_distribution.dat
ncolumns=`cat ${rdf_file} | awk '{ if(NR==1) print NF }'`

i=1
while [ $i -le ${ncolumns} ]
do
  awk '{print $'$i'}' ${rdf_file} | blocked_error 
  i=`expr $i + 1`
done

