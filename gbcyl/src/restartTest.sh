#!/bin/bash


# Makes various tests to the restart feature of gbcyl.


# 0. Preparations
mkdir -p testtmp
cp ../data/r9-n1290-cr.gb testtmp/ 
cp ../data/gbcyl.in testtmp/gbcyl.in.tmp
cp gbcyl testtmp/
cd testtmp
test -f simdata.out && rm -f simdata.out


# 1. Test that restarting does not change anything.

# 1.1. Equilibrate 0 sweeps. Make production for 0 sweeps.
cat gbcyl.in.tmp| 
awk '
  /Nprod/{$2 = "0"} 
  /Nrelax/{$2 = "0"} 
  {print}' > gbcyl.in
./gbcyl

# 1.2. Save simdata.out as simdata.out.old
mv simdata.out simdata.out.old

# 1.3. Save restart.gbcyl as restart.gbcyl.old
cp restart.gbcyl restart.gbcyl.old

# 1.4. Restart.
./gbcyl restart 

# 1.5. Compare simdata.out with simdata.out.old
diff -q simdata.out simdata.out.old 
DIFF_SIMDATA=$?

# 1.6. Compare restart.gbcyl with restart.gbcyl.old
diff -q restart.gbcyl restart.gbcyl.old
DIFF_RESTARTFILE=$?

# 1.7. Clean up. 
rm *.old
rm simdata.out

# 1.8. Exit if failed. 
if [ "${DIFF_SIMDATA}" -ne 0 ]
then
  echo "Restart changed file simdata.out." 
  exit 1
fi
if [ "${DIFF_RESTARTFILE}" -ne 0 ]
then
  echo "Restart changed file restart.gbcyl"
  exit 1
fi



# 2. Test that restarting results in the same configurations of molecules

# 2.1. Edit restart.gbcyl: production 10 sweeps, current sweep = 0.
cp restart.gbcyl restart.gbcyl.old
cat restart.gbcyl.old| 
awk '
  /N_PRODUCTION/{$2 = "10,"}
  /N_EQUILIBRATION/{$2 = "0,"}
  /I_SWEEP/{$2 = "0,"} 
  {print}' > restart.gbcyl

# 2.2. Restart simulation.
./gbcyl restart

# 2.3. Save simdata.out as simdata.out.old
mv simdata.out simdata.out.old

# 2.4. Save restart.gbcyl as restart.gbcyl.old
cp restart.gbcyl restart.gbcyl.old

# 2.5. Edit gbcyl.in: Production sweeps = 10 
cat gbcyl.in.tmp| awk '
  /Nprod/{$2 = "10"}
  /Nrelax/{$2 = "0"}
  {print}' > gbcyl.in

# 2.6. Run simulation.
./gbcyl

# 2.7. Compare simdata.out with simdata.out.old
diff -q simdata.out simdata.out.old
DIFF_SIMDATA=$?

# 2.8. Compare restart.gbcyl with restart.gbcyl.old
diff -q restart.gbcyl restart.gbcyl.old
DIFF_RESTARTFILE=$?

# 2.9. Clean up.
rm -f *.old
rm -f simdata.out 

# 2.10. Exit if failed.
if [ "${DIFF_SIMDATA}" -ne 0 ]
then
  echo "simdata.out different when initiating with restart file."
  exit 1
fi
if [ "${DIFF_RESTARTFILE}" -ne 0 ]
then
  echo "restart.gbcyl different when initiating with restart file." 
  exit 1
fi


# 3. Test that continuing a run works the same way as running at once.

# 3.1. Run equilibration 5 sweeps, production 10 sweeps.
cat gbcyl.in.tmp|
awk '
  /Nprod/{$2 = "10"}
  /Nrelax/{$2 = "5"}
  {print}' > gbcyl.in
./gbcyl

# 3.2. Edit restart.gbcyl: Production sweeps = 15.
cp restart.gbcyl restart.gbcyl.old
cat restart.gbcyl.old|
awk '
  /N_PRODUCTION/{$2 = "15,"}
  {print}' > restart.gbcyl

# 3.3. Restart.
./gbcyl restart

# 3.4. Save simdata.out as simdata.out.old.
mv simdata.out simdata.out.old

# 3.5. Save restart.gbcyl as restart.gbcyl.old.
cp restart.gbcyl restart.gbcyl.old

# 3.6. Run equilibration 5 sweeps, production 15 sweeps.
cat gbcyl.in.tmp|
awk '
  /Nrelax/{$2 = "5"}
  /Nprod/{$2 = "15"}
  {print}' > gbcyl.in
./gbcyl

# 3.7. Compare simdata.out with simdata.out.old.
diff -q simdata.out simdata.out.old
DIFF_SIMDATA=$?

# 3.8. Compare restart.gbcyl with restart.gbcyl.old.
diff -q restart.gbcyl restart.gbcyl.old
DIFF_RESTARTFILE=$?

# 3.9. Clean up. 
rm -f simdata.out
rm -f *.old

# 3.10. Exit if failed.
if [ "${DIFF_SIMDATA}" -ne 0 ]
then
  echo "Running at once is different from continuing in simdatas."
  exit 1
fi
if [ "${DIFF_RESTARTFILE}" -ne 0 ]
then 
  echo "Running at once is different from continuing in restart files."
  exit 1
fi


# 4. Clean up. 
# :NOTE: This won't work if a preceding test fails. 
cd ..
rm -rf testtmp