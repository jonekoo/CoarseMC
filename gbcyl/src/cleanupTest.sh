#!/bin/sh

if [ ${OMPI_COMM_WORLD_RANK} == 0 ]
then 
  rm fin-kgb.*.* ptgbcyl.in r9-n1290-cr.gb mt-state.*
fi
