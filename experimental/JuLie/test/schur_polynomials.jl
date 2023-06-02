import Oscar.JuLie: schur_polynomial_cbf

@testset "Schur polynomials" begin
	S,x = polynomial_ring(ZZ, :x => 1:3)

	# schur_polynomial(R, λ)
	@test schur_polynomial(S, Partition([1])) == x[1]
	@test schur_polynomial(S, Partition([1,1])) == x[1]*x[2]
	@test schur_polynomial(S, Partition([2])) == x[1]^2
	@test schur_polynomial(S, Partition([1,1,1])) == x[1]*x[2]*x[3]
	@test schur_polynomial(S, Partition([2,1])) == x[1]^2*x[2] + x[1]*x[2]^2
	@test schur_polynomial(S, Partition([3])) == x[1]^3
	@test schur_polynomial(S, Partition([])) == 1

	# schur_polynomial(R, λ, n)
	@test schur_polynomial(S, Partition([1]), 3) == x[1] + x[2] + x[3]
	@test schur_polynomial(S, Partition([1,1]), 3) == x[1]*x[2] + x[1]*x[3] + x[2]*x[3]
	@test schur_polynomial(S, Partition([2]), 2) == x[1]^2 + x[2]^2 + x[1]*x[2]
	@test schur_polynomial(S, Partition([2]), 3) == x[1]^2 + x[2]^2 + x[1]*x[2] + x[1]*x[3] + x[2]*x[3] + x[3]^2
	@test schur_polynomial(S, Partition([1,1,1]), 3) == x[1]*x[2]*x[3]
	@test schur_polynomial(S, Partition([2,1]), 1) == 0
	@test schur_polynomial(S, Partition([2,1]), 3) == 2*x[1]*x[2]*x[3] + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[3] + x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2
	@test schur_polynomial(S, Partition([3,1]), 2) == x[1]^3*x[2] + x[1]^2*x[2]^2 + x[1]*x[2]^3
	@test schur_polynomial(S, Partition([3]), 2) == x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^3 + x[2]^3
	@test schur_polynomial(S, Partition([3]), 3) == x[1]*x[2]*x[3] + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[3] + x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2 + x[1]^3 + x[2]^3 + x[3]^3

	@test schur_polynomial(S, Partition([141]), 1) == x[1]^141 #this calls Cauchy's algorithm
	@test schur_polynomial(S, Partition([200]), 1) == x[1]^200 #this calls Cauchy's algorithm

	@test schur_polynomial(S, Partition([]), 1) == 1
	@test schur_polynomial(S, Partition([]), 0) == 1
	@test schur_polynomial(S, Partition([3,2,1]), 0) == 0
	@test_throws ArgumentError schur_polynomial(Partition([3,2,1]), -1)

	#schur_polynomial(λ)
	@test schur_polynomial(Partition([])) == 1
	@test schur_polynomial(Partition([1]))(1) == 1

	#schur_polynomial(λ, n)
	@test schur_polynomial(Partition([]), 1) == 1
	@test schur_polynomial(Partition([]), 0) == 1
	@test schur_polynomial(Partition([3,2,1]), 0) == 0
	@test schur_polynomial(Partition([1]), 10)(1,1,1,1,1,1,1,1,1,1) == 10
	@test_throws ArgumentError schur_polynomial(Partition([3,2,1]), -1)

	#Two examples from Wikipedia
	schur_polynomial(S, Partition([2,1,1]), 3) == x[1]*x[2]*x[3]*(x[1]+x[2]+x[3])
	schur_polynomial(S, Partition([2,2]), 3) == x[1]^2*x[2]^2 + x[1]^2*x[3]^2 + x[2]^2*x[3]^2 + x[1]^2*x[2]*x[3] + x[1]*x[2]^2*x[3] + x[1]*x[2]*x[3]^2

end
