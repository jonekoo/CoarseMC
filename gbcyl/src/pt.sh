#!/bin/sh

#--------------------------------------
# copy files to rundir
##RUNDIR=${HOME}/test/rundir
RUNDIR=$1
mkdir -p ${RUNDIR}
cp fin-kgb.${OMPI_COMM_WORLD_RANK}.pms \
fin-kgb.${OMPI_COMM_WORLD_RANK}.mol ${RUNDIR}/

#---------------
# run simulation
cd ${RUNDIR}
ptgbcyl2.3


