#!/bin/sh

INSTALLDIR=${HOME}

PARAMETERFILE=${INSTALLDIR}/share/coarsemc.in
GEOMETRYFILE=${INSTALLDIR}/share/r9-n1290-cr.gb

i=0
while [ $i -lt $1 ]
do
  cp ${GEOMETRYFILE} fin-kgb.${i}.mol

  #------------------------------------
  # calculate temperature for this task
  tempmin=$2
  tempstep=$3
  temp=$(echo "${tempmin} + ${i} * ${tempstep}" |bc -l)

  #------------------------------
  # set temperature for this task
  cat ${PARAMETERFILE} | \
  awk -v temp=${temp} '$1 ~ /temperature/ {$2=temp}; {print}' > \
  fin-kgb.${i}.pms

  #-----------------------------
  # fix the geometry file format
  ##sed -i 's/ $R/$R/g' fin-kgb.${i}.mol
  i=`expr $i + 1`
done
