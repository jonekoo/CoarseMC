#!/bin/bash

# Makes various tests to the restart feature of ptgbcyl.

sh ${srcdir}/initTest.sh

id=${OMPI_COMM_WORLD_RANK}

# 1. 
if [ "${id}" = "0" ]
then
  echo "Test that restarting does not change anything."
fi

# 1.1. Equilibrate 0 sweeps. Make production for 0 sweeps.
cat ptgbcyl.in| \
awk '
  /n_equilibration_sweeps/{$2 = "0"} 
  /n_production_sweeps/{$2 = "0"} 
  {print}' > fin-kgb.${id}.pms
cp r9-n1290-cr.gb fin-kgb.${id}.mol

./ptgbcyl >> ptgbcyl.${id}.out

# 1.2. Copy current molecule and parameter files for later comparison
cp fin-kgb.${id}.mol fin-kgb.${id}.mol.old
cp fin-kgb.${id}.pms fin-kgb.${id}.pms.old

# 1.3. Restart.
./ptgbcyl >> ptgbcyl.${id}.out

# 1.4. Compare molecule files
diff -q fin-kgb.${id}.mol fin-kgb.${id}.mol.old 
DIFF_SIMDATA=$?

# 1.5. Compare parameter files
diff -q fin-kgb.${id}.pms fin-kgb.${id}.pms.old
DIFF_RESTARTFILE=$?

# 1.7. Exit if failed. 
if [ "${DIFF_SIMDATA}" -ne 0 ]
then
  echo "Restart changed file fin-kgb.${id}.mol." 
  sh ${srcdir}/cleanupTest.sh
  exit 1
fi
if [ "${DIFF_RESTARTFILE}" -ne 0 ]
then
  echo "Restart changed file fin-kgb.${id}.pms"
  sh ${srcdir}/cleanupTest.sh
  exit 1
fi


sh ${srcdir}/cleanupTest.sh