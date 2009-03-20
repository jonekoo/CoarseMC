bin_PROGRAMS = gbcyl energy_calculator
gbcyl_SOURCES = nrtype.f90 nrutil.f90 nr.f90 attwall.f90 repwall.f90\
  cylinder.f90 gayberne.f90 mtfort90.f90 utils.f90 particle.f90\
  particlewall.f90\
  verlet.f90 energy.f90 horner.f90 io.f90 mc_sweep.f90 \
  rd_s.f90 rd_v.f90 rf_s.f90 rf_v.f90 mc_engine.f90\
  gbcyl.f90

energy_calculator_SOURCES = nrtype.f90 \
                 nrutil.f90 \
                 nr.f90 \
                 utils.f90 \
                 mtfort90.f90 \
                 gayberne.f90 \
                 particle.f90 \
                 state_reader.f90 \
                 io.f90 \
                 horner.f90 \
                 attwall.f90 \
                 repwall.f90 \
                 cylinder.f90 \
                 particlewall.f90 \
                 verlet.f90 \
                 energy.f90 \
                 energy_calculator.f90 \
                 rd_s.f90 \
                 rd_v.f90 \
                 rf_s.f90 \
                 rf_v.f90 

TESTS = TestRunner
dist_check_SCRIPTS = restartTest.sh
check_PROGRAMS = TestRunner
TestRunner_SOURCES = nrtype.f90 gayberne.f90 gayberne_fun.f90 TestRunner.f90
MODEXT = @MODEXT@


mostlyclean-local :
	-rm -f *.$(MODEXT)
