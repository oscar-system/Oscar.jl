@testset "combinatorics/q_analogs.jl" begin
  R,q = laurent_polynomial_ring(ZZ, "q")

	# Test if everything specializes to usual stuff at q=1
	@test q_integer(0,1) == 0
	@test q_integer(1,1) == 1
	@test q_integer(2,1) == 2
	@test q_integer(3,1) == 3
	@test q_integer(-1,1) == -1
	@test q_integer(-2,1) == -2
	@test q_integer(-3,1) == -3

	@test q_factorial(0,1) == factorial(0)
	@test q_factorial(1,1) == factorial(1)
	@test q_factorial(2,1) == factorial(2)
	@test q_factorial(3,1) == factorial(3)

	@test q_binomial(0,0,1) == binomial(0,0)
	@test q_binomial(0,1,1) == binomial(0,1)
	@test q_binomial(0,2,1) == binomial(0,2)
	@test q_binomial(0,3,1) == binomial(0,3)
	@test q_binomial(1,1,1) == binomial(1,1)
	@test q_binomial(2,1,1) == binomial(2,1)
	@test q_binomial(3,1,1) == binomial(3,1)
	@test q_binomial(2,2,1) == binomial(2,2)
	@test q_binomial(3,2,1) == binomial(3,2)
	@test q_binomial(-3,0,1) == binomial(-3,0)
	@test q_binomial(-3,1,1) == binomial(-3,1)
	@test q_binomial(-3,2,1) == binomial(-3,2)
	@test q_binomial(2,3,1) == binomial(2,3)

	# Examples from Conrad (2000)
	@test q_integer(0) == 0
	@test q_integer(1) == 1
	@test q_integer(2) == 1+q
	@test q_integer(-1) == -q^(-1)

	@test q_factorial(0) == 1
	@test q_factorial(1) == 1
	@test q_factorial(2) == 1+q
	@test q_factorial(3) == 1+2*q + 2*q^2 + q^3

	@test q_binomial(0,0) == 1
	@test q_binomial(1,0) == 1
	@test q_binomial(2,0) == 1

	@test q_binomial(0,1) == q_integer(0)
	@test q_binomial(1,1) == q_integer(1)
	@test q_binomial(2,1) == q_integer(2)

	@test q_binomial(4,2) == (1+q^2)*(1+q+q^2)

	K,i = cyclotomic_field(4);
	@test q_binomial(4,2,i) == 0

	@test q_binomial(19,5,-1) == 36
	@test q_binomial(17,10,i) == 0
	@test q_binomial(-5,6,i) == -2*i

	# From https://en.wikipedia.org/wiki/Gaussian_binomial_coefficient
	@test q_binomial(0,0) == 1
	@test q_binomial(1,0) == 1
	@test q_binomial(1,1) == 1
	@test q_binomial(2,1) == 1+q
	@test q_binomial(3,1) == 1+q+q^2
	@test q_binomial(3,2) == 1+q+q^2
	@test q_binomial(4,2) == (1+q^2)*(1+q+q^2)
	@test q_binomial(6,3) == (1+q^2)*(1+q^3)*(1+q+q^2+q^3+q^4)

end
