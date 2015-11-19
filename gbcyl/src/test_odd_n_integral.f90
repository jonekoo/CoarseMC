program test
use odd_n_integral
use iso_fortran_env, only: dp => REAL64
implicit none

write(*, *) "0! = ", factorial(0)
write(*, *) "1! = ", factorial(1)
write(*, *) "3! = ", factorial(3)
write(*, *) ""
write(*, *) "(1/2)_0 = ", pochhammer(0.5_dp, 0)
write(*, *) "(0)_3 = ", pochhammer(0._dp, 3)
write(*, *) "(0)_0 = ", pochhammer(0._dp, 0)
write(*, *) "(1/2)_3 = ", pochhammer(0.5_dp, 3)
write(*, *) "(2)_3 = ", pochhammer(2._dp, 3)
write(*, *) ""
write(*, *) "binomial(0, 0) = ", binomial(0, 0)
write(*, *) "binomial(1, 0) = ", binomial(1, 0)
write(*, *) "binomial(0, 1) = ", binomial(0, 1)
write(*, *) "binomial(5, 2) = ", binomial(5, 2)
write(*, *) ""
write(*, *) "I_8^0 = ", angle_integral(8, 0, 1._dp)
write(*, *) "I_36^2 = ", angle_integral(36, 2, 0.95_dp)
write(*, *) ""
write(*, *) "angle_integral2(8, 0, 0.9_dp) = ", angle_integral2(8, 0, 0.9_dp)
write(*, *) "angle_integral2(7, 0, 0.9_dp) = ", angle_integral2(7, 0, 0.9_dp)
write(*, *) "angle_integral(7, 0, 0.9_dp) = ", angle_integral(7, 0, 0.9_dp)
write(*, *) ""
write(*, *) "innersum(8,0,0,1) = ", innersum(8, 0, 0, 1._dp)
write(*, *) "innersum(8,2,2,1) = ", innersum(8, 2, 2, 1._dp)
write(*, *) "innersum(8,2,4,0.9) = ", innersum(8, 2, 4, 0.9_dp)
write(*, *) "innersum(7,0,0,0.9) = ", innersum(7, 0, 0, 0.9_dp)
write(*, *) ""
write(*, *) "integral(n=11, m=0, k=0.9, R=9.0) = ", integral(11, 0, 0.9_dp, 9._dp)
write(*, *) "integral(n=37, m=0, k=0.9, R=9.0) = ", integral(37, 0, 0.9_dp, 9._dp)
write(*, *) "integral(n=11, m=2, k=0.9, R=9.0) = ", integral(11, 2, 0.9_dp, 9._dp)
write(*, *) "integral(n=37, m=2, k=0.9, R=9.0) = ", integral(37, 2, 0.95_dp, 9._dp)
end program
