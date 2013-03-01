program testrunner
use ftnunit
use gayberne_test, only: gb_test_all => test_all 
use particle_test, only: particle_test_all => test_all
use class_poly_box_test, only: box_test_all => test_all
use class_factory_test, only: factory_test_all => test_all
use utils_test, only: utils_test_all => test_all
implicit none

call runtests_init
call runtests(gb_test_all)
call runtests(particle_test_all)
call runtests(box_test_all)
call runtests(factory_test_all)
call runtests(utils_test_all)
call runtests_final

end program
