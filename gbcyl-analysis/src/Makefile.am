bin_PROGRAMS = orientation_parameter \
               cylindrical_orientation \
               tau1 \
               radial_distribution \
               density \
               psi6 \
               povmaker \
               particlehistogram 

dist_bin_SCRIPTS = rdf_average.sh \
                   loop-script-over-args.sh \
                   remove-non-numeric-lines.sh \
                   density_cyl.awk \
                   show-non-numeric-lines.sh \
                   select.awk \
                   periodic_select.awk \
                   volume_cyl.awk \
                   find_closest.awk \
                   dots2commas.sh \
                   block-by-column.sh

orientation_parameter_SOURCES = nrtype.f90 utils.f90 particle.f90 \
state_reader.f90 orientational_ordering.f90 orientation_parameter_wrapper.f90 \
orientation_parameter.f90

cylindrical_orientation_SOURCES = nrtype.f90 utils.f90 particle.f90 \
state_reader.f90 orientational_ordering.f90 \
cylindrical_orientation.f90 cylindrical_orientation_wrapper.f90 \
cylindrical_orientation_tensor.f90

tau1_SOURCES = nrtype.f90 nrutil.f90 mnbrak.f90 brent.f90 nr.f90 \
orientational_ordering.f90 tau1.f90 tau1_module.f90 state_reader.f90 \
tau1_program.f90 

radial_distribution_SOURCES = nrtype.f90 \
                              utils.f90 \
                              particle.f90 \
                              state_reader.f90 \
                              radial_distribution_wrapper.f90 \
                              radial_distribution.f90 \
                              gr3d.f90

density_SOURCES = nrtype.f90 \
                  utils.f90 \
                  particle.f90 \
                  state_reader.f90 \
                  density.f90

psi6_SOURCES = nrtype.f90 \
               utils.f90 \
               particle.f90 \
               state_reader.f90 \
               cylinder.f90 \
               geometry.f90 \
               psi6_module.f90 \
               orientational_ordering.f90 \
               psi6.f90

povmaker_SOURCES = nrtype.f90 \
                   utils.f90 \
                   particle.f90 \
                   io.f90 \
                   state_reader.f90 \
                   povize.f90 
                   
particlehistogram_SOURCES = nrtype.f90 \
                            utils.f90 \
                            particle.f90 \
                            state_reader.f90 \
                            histogram.f90 \
                            orientational_ordering.f90 \
                            particlehistogram.f90 \
                            make_bin.f90 \
                            bin_indices.f90 \
                            rho.f90


MODEXT = @MODEXT@

mostlyclean-local :
	-rm -f *.$(MODEXT)

