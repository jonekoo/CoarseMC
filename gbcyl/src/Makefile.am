bin_PROGRAMS = gbcyl
gbcyl_SOURCES = nrtype.f90 nrutil.f90 nr.f90 attwall.f90 repwall.f90\
  cylinder.f90 gayberne.f90 mtfort90.f90 utils.f90 particle.f90\
  particlewall.f90\
  verlet.f90 energy.f90 horner.f90 io.f90 mc_sweep.f90 \
  nr.f90 nrtype.f90 rd_s.f90 rd_v.f90 rf_s.f90 rf_v.f90 mc_engine.f90\
  gbcyl.f90 npt_sweep.f90

MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)
