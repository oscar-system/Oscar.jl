import .JuLie: schur_polynomial_cbf

@testset "Schur polynomials" begin
	x = [string("x",string(i)) for i=1:3]
	S,x = PolynomialRing(ZZ, x)

	# schur_polynomial(λ, R)
	@test schur_polynomial(Partition([1]), S) == x[1]
	@test schur_polynomial(Partition([1,1]), S) == x[1]*x[2]
	@test schur_polynomial(Partition([2]), S) == x[1]^2 + x[2]^2 + x[1]*x[2]
	@test schur_polynomial(Partition([1,1,1]), S) == x[1]*x[2]*x[3]
	@test schur_polynomial(Partition([2,1]), S) == 2*x[1]*x[2]*x[3] + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[3] + x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2
	@test schur_polynomial(Partition([3]), S) == x[1]*x[2]*x[3] + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[3] + x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2 + x[1]^3 + x[2]^3 + x[3]^3
	@test schur_polynomial(Partition([]), S) == 1

	# schur_polynomial(λ, R, n)
	@test schur_polynomial(Partition([1]), S, 3) == x[1] + x[2] + x[3]
	@test schur_polynomial(Partition([1,1]), S, 3) == x[1]*x[2] + x[1]*x[3] + x[2]*x[3]
	@test schur_polynomial(Partition([2]), S, 3) == x[1]^2 + x[2]^2 + x[1]*x[2] + x[1]*x[3] + x[2]*x[3] + x[3]^2
	@test schur_polynomial(Partition([1,1,1]), S, 3) == x[1]*x[2]*x[3]
	@test schur_polynomial(Partition([1,1,1]), S, 4) == x[1]*x[2]*x[3]
	@test schur_polynomial(Partition([2,1]), S, 1) == 0
	@test schur_polynomial(Partition([3,1]), S, 2) == x[1]^3*x[2] + x[1]^2*x[2]^2 + x[1]*x[2]^3
	@test schur_polynomial(Partition([3]), S, 2) == x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^3 + x[2]^3

	@test schur_polynomial(Partition([141]), S, 1) == x[1]^141 #this calls Cauchy's algorithm
	@test schur_polynomial(Partition([200]), S, 1) == x[1]^200 #this calls Cauchy's algorithm

	@test schur_polynomial(Partition([]), S, 1) == 1
	@test schur_polynomial(Partition([]), S, 0) == 1
	@test schur_polynomial(Partition([3,2,1]), S, 0) == 0
	@test_throws ArgumentError schur_polynomial(Partition([3,2,1]), -1)

	# schur_polynomial(λ, x)
	@test schur_polynomial(Partition([]), [x[1]]) == 1
	@test schur_polynomial(Partition([]), x[1:0]) == 1
	@test schur_polynomial(Partition([1]), x[1:0]) == 0
	@test schur_polynomial(Partition([141]), [x[1]]) == x[1]^141 #this calls Cauchy's algorithm
	@test schur_polynomial(Partition([200]), [x[2]]) == x[2]^200 #this calls Cauchy's algorithm
	@test schur_polynomial(Partition([1]), PolynomialRing(ZZ, [string("x",string(i)) for i=1:10])[2])(1,1,1,1,1,1,1,1,1,1) == 10 #this calls Cauchy's algorithm
	@test schur_polynomial_cbf(Partition([3,2,1]), [x[1],x[2],x[3]]) == schur_polynomial(Partition([3,2,1]), [x[1],x[2],x[3]])

	check = true
	for n = 1:5
		x = [string("x",string(i)) for i=1:n]
		S,x = PolynomialRing(ZZ, x)
		for p in partitions(n)
			if schur_polynomial(p, x) != schur_polynomial(p, S)
				check = false
			end
		end
	end
	@test check == true

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
	schur_polynomial(Partition([2,1,1]), [x[1],x[2],x[3]]) == x[1]*x[2]*x[3]*(x[1]+x[2]+x[3])
	schur_polynomial(Partition([2,2]), [x[1],x[2],x[3]]) == x[1]^2*x[2]^2 + x[1]^2*x[3]^2 + x[2]^2*x[3]^2 + x[1]^2*x[2]*x[3] + x[1]*x[2]^2*x[3] + x[1]*x[2]*x[3]^2

end