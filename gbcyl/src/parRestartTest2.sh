#!/bin/sh


sh ${srcdir}/initTest.sh

#cd testtmp
#cd ${srcdir}

id=${OMPI_COMM_WORLD_RANK}

# 2. 
if [ ${id} = 0 ] 
then
  echo "Test that restarting results in the same configuration of molecules"
  echo "First run 0 + 1 sweeps in two simulations, then run 1 sweep in one"
  echo "simulation. Test is passed only if resulting parameter and molecule"
  echo "files are exactly the same."
fi

cat ptgbcyl.in| \
awk '
  /n_production_sweeps/{$2 = "0"}
  /n_equilibration_sweeps/{$2 = "0"}
  /production_period/{$2 = "1"}
  /i_sweep/{$2 = "0"} 
  {print}' > fin-kgb.${id}.pms
cp r9-n1290-cr.gb fin-kgb.${id}.mol

# 2.2. Start simulation.
./ptgbcyl >> ptgbcyl.${id}.out

cp fin-kgb.${id}.pms fin-kgb.${id}.pms.tmp
cat fin-kgb.${id}.pms.tmp | \
awk '
  /n_production_sweeps/{$2 = "1"}
  {print}' > fin-kgb.${id}.pms

# 2.3 Restart simulation
./ptgbcyl >> ptgbcyl.${id}.out

##nconfs=`grep 'R' fin-kgb.${id}.mol | wc -l`
##echo "molecule file has ${nconfs} configurations"

# 2.3. Save copy of current molecule file before restarting
mv fin-kgb.${id}.mol fin-kgb.${id}.mol.old

# 2.4. Save copy of current parameter file before restarting
mv fin-kgb.${id}.pms fin-kgb.${id}.pms.old

# 2.5 Rename rng state file so that it won't be used in the next run.
mv mt-state.${id} mt-state.${id}.old

# 2.6. Edit ptgbcyl.in: Production sweeps = 1 
cat ptgbcyl.in| awk '
  /n_production_sweeps/{$2 = "1"}
  /n_equilibration_sweeps/{$2 = "0"}
  /production_period/{$2 = "1"}
  {print}' > fin-kgb.${id}.pms
cp r9-n1290-cr.gb fin-kgb.${id}.mol

# 2.7. Run simulation.
./ptgbcyl >> ptgbcyl.${id}.out

# 2.8. Compare fin-kgb.${id}.mol with fin-kgb.${id}.mol.old
diff -q fin-kgb.${id}.mol fin-kgb.${id}.mol.old
DIFF_SIMDATA=$?

# 2.9. Compare fin-kgb.${id}.pms with fin-kgb.${id}.pms.old
diff -q fin-kgb.${id}.pms fin-kgb.${id}.pms.old
DIFF_RESTARTFILE=$?

# 2.10. Exit if failed.
if [ "${DIFF_SIMDATA}" -ne 0 ]
then
  echo "fin-kgb.${id}.mol different when initiating with restart file."
  sh ${srcdir}/cleanupTest.sh
  exit 1
fi
if [ "${DIFF_RESTARTFILE}" -ne 0 ]
then
  echo "fin-kgb.${id}.pms different when initiating with restart file." 
  sh ${srcdir}/cleanupTest.sh
  exit 1
fi

sh ${srcdir}/cleanupTest.sh
