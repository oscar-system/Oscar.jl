using Random

RNG = Random.MersenneTwister(42)

"""
	randpoly(R::Oscar.CRing,coeffs=0:9,max_exp=4,max_terms=8)
> Return a random Polynomial from the Polynomial Ring `R` with coefficients in `coeffs`
> with exponents between `0` and `max_exp` und between `0` and `max_terms` terms
"""
function randpoly(R::Oscar.CRing,coeffs=0:9,max_exp=4,max_terms=8)
	n = nvars(R)
	K = base_ring(R)
	E = [[Random.rand(RNG,0:max_exp) for i=1:n] for j=1:max_terms]
	C = [K(Random.rand(RNG,coeffs)) for i=1:max_terms]
	M = MPolyBuildCtx(R)
	for i=1:max_terms
		push_term!(M,C[i],E[i])
	end
	return finish(M)
end

"""
	array_to_matrix(A::Array,R::CRing = parent(A[1,1]))
Return `A` as an AbstractAlgebra Matrix
"""
function array_to_matrix(A::Array,R::AbstractAlgebra.Ring = parent(A[1,1]))
	Mat = AbstractAlgebra.MatrixSpace(R,size(A)...)
	return Mat(R.(A))
end

@testset "Test intersection" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

  A1 = R[x y;
        2*x^2 3*y^2]
  A2 = R[x^3 x^2*y;
        (2*x^2+x) (2*y^3+y)]
  B = R[4*x*y^3 (2*x+y)]

  M1 = Oscar.SubQuo(A1, B)
  M2 = Oscar.SubQuo(A2, B)

  A1 = R[x R(1); y x^2]
  A2 = R[x y]
  B1 = R[x^2 y^2]
  res = R[(x*y^2) (y^3-x*y+y^2); (-x*y) (-x^2*y^2+x*y^3-x^2*y+y^3-x*y)]

  @test M1 == intersect(M1, M1)[1]
  @test Oscar.SubQuo(res, B1) == intersect(Oscar.SubQuo(A1, B1), Oscar.SubQuo(A2, B1))[1] #macaulay2

  A1 = R[x y; x^2 y^2]
  A2 = R[x+x^2 y+y^2]
  @test Oscar.SubQuo(A2, B1) == intersect(Oscar.SubQuo(A1,B1), Oscar.SubQuo(A2,B1))[1]
end

@testset "Test presentation" begin

	# over Integers
	#=R, (x,y,z) = PolynomialRing(ZZ, ["x", "y", "z"])
	generator_matrices = [R[x x^2*y; y y^2*x^2], R[x x^2*y; y y^2*x^2], R[x y; x^2*y y^2*x], R[x+R(1) x^10; x^2+y^2 y^4-x^4], R[42*x*y 7*x^2; 6*x 9*y^2]]
	relation_matrices  = [R[x^3 y^4], R[x^3 y^4; x^2*y x*y^2], R[x*y^2 x; y^3 x*y^3], R[x x*y], R[3*x*y 7*x*y; 42 7]]
	true_pres_matrices = [R[x^5*y-y^4 -x^5+x*y^3], R[-x^2*y-x*y-y x^2+x; x*y^4-x^3*y+y^4-x^2*y -x*y^3+x^3; -y^4+x^2*y x*y^3-x^3; x^2*y^3+x*y^3+y^3 -x^2*y^2-x*y^2], R[-x*y R(1); -x^2*y^5+x*y^3 R(0)], R[x^5-x*y^4+x^3*y+x*y^3 x^11-x^2*y-x*y; -x^5*y^7+x*y^11+2*x^5*y^6-x^3*y^8-3*x*y^10+2*x^5*y^5+2*x^3*y^7-3*x^5*y^4+2*x^3*y^6+5*x*y^8+2*x^5*y^3-3*x^3*y^5-5*x*y^7+2*x^3*y^4+2*x*y^6 -x^11*y^7+2*x^11*y^6+2*x^11*y^5-3*x^11*y^4+2*x^11*y^3+x^2*y^8-2*x^2*y^7+x*y^8-2*x^2*y^6-2*x*y^7+3*x^2*y^5-2*x*y^6-2*x^2*y^4+3*x*y^5-2*x*y^4], R[-13*y R(0); -377 -2639*x; -39 -273*x; -13*y-39 -273*x; -13 -91*x; y^2-42*x -294*x^2+21*x*y; 9*y^2-x+26*y+78 -7*x^2+189*x*y+546*x; -y^2+3*x 21*x^2-21*x*y]]
	for (A,B,true_pres_mat) in zip(generator_matrices, relation_matrices, true_pres_matrices)
    println(A,B,true_pres_mat)
		SQ = Oscar.SubQuo(A,B)
		pres_mat = matrix(Oscar.present_as_cokernel(SQ).quo)
		@test Oscar.cokernel(pres_mat) == Oscar.cokernel(true_pres_mat)

		pres_SQ, i = Oscar.present_as_cokernel(SQ, :both)
		p = i.inverse_isomorphism
		@test Oscar.iswelldefined(i)
		@test Oscar.iswelldefined(p)
		@test Oscar.isbijective(i)
		@test Oscar.isbijective(p)
	end=#




	# over Rationals
	R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	generator_matrices = [R[x x^2*y; y y^2*x^2], R[x x^2*y; y y^2*x^2], R[x y; x^2*y y^2*x], R[x+R(1) x^10; x^2+y^2 y^4-x^4], R[x+R(1) x^10; x^2+y^2 y^4-x^4]]
	relation_matrices  = [R[x^3 y^4], R[x^3 y^4; x^2*y x*y^2], R[x*y^2 x; y^3 x*y^3], R[x+y x+y; x*y x*y], R[x+y x+y; y x; x y]]
	true_pres_matrices = [R[x^5*y-y^4 -x^5+x*y^3], R[-x^2*y-x*y-y x^2+x; -x^2*y^4+x^4*y x^2*y^3-x^4], R[-x*y R(1); -x^2*y^5+x*y^3 R(0)], R[-x^4+y^4-x^2-y^2 -x^10+x+R(1)], R[R(0) -x^2+y^2; -2*x^2 -x^9*y+x+1; -2*x*y^2 -x^8*y^3+y^2+x; -2*x*y^2+2*y^2 -x^8*y^3+x^7*y^3+y^2-1]]
	for (A,B,true_pres_mat) in zip(generator_matrices, relation_matrices, true_pres_matrices)
		SQ = Oscar.SubQuo(A,B)
		pres_mat = matrix(Oscar.present_as_cokernel(SQ).quo)
		@test Oscar.cokernel(pres_mat) == Oscar.cokernel(true_pres_mat)

		pres_SQ, i= Oscar.present_as_cokernel(SQ, :both)
    p = i.inverse_isomorphism
		@test Oscar.iswelldefined(i)
		@test Oscar.iswelldefined(p)
		@test Oscar.isbijective(i)
		@test Oscar.isbijective(p)
	end
end


#=@testset "Test kernel" begin
	# This test doesn't terminate due to a bug in syz

	# over Integers
	R, (x,y,z) = PolynomialRing(ZZ, ["x", "y", "z"])
	matrices = [R[x^2+x y^2+y; x^2+y y^2; x y], R[5*x^5+x*y^2 4*x*y+y^2+R(1); 4*x^2*y 2*x^2+3*y^2-R(5)]]
	kernels = [R[x^2-x*y+y -x^2+x*y -x^2+x*y-y^2-y], R[R(0) R(0)]]
	for (A,Ker) in zip(matrices, kernels)
		K,emb = Oscar.kernel(Oscar.matrix_to_map(A))
    @test K == Oscar.image(emb)[1]
		#@test image(K) == image(Ker)
    @test Oscar.image(emb)[1] == Oscar.image(Oscar.matrix_to_map(Ker))[1]
	end
	#=for k=1:3
		A = array_to_matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
    A = Oscar.matrix_to_map(A)
		K,emb = Oscar.kernel(A)
		@test Oscar.iszero(emb*A)
	end=#

	# over Rationals
	R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	matrices = [R[x^2+x y^2+y; x^2+y y^2; x y], R[5*x^5+x*y^2 4*x*y+y^2+R(1); 4*x^2*y 2*x^2+3*y^2-R(5)]]
	kernels = [R[x^2-x*y+y -x^2+x*y -x^2+x*y-y^2-y], R[R(0) R(0)]]
	for (A,Ker) in zip(matrices, kernels)
		K,emb = Oscar.kernel(Oscar.matrix_to_map(A))
    @test K == Oscar.image(emb)[1]
    @test Oscar.image(emb)[1] == Oscar.image(Oscar.matrix_to_map(Ker))[1]
	end
	for k=1:3
		A = array_to_matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
    display(A)
    A = Oscar.matrix_to_map(A)
		K,emb = Oscar.kernel(A)
    @test Oscar.image(emb)[1] == K
		@test Oscar.iszero(emb*A)
	end
end=#

@testset "Test iszero(SubQuo)" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	A = R[x^2+2*x*y y^2*x-2*x^2*y;-y x*y]
	B = R[x^2 y^2*x;-y x*y]
	@test Oscar.iszero(Oscar.SubQuo(A,B))
	for k=1:3
		A = array_to_matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
		B = array_to_matrix([randpoly(R,0:15,2,2) for i=1:2,j=1:2])
		@test Oscar.iszero(Oscar.SubQuo(A,A))
		@test !Oscar.iszero(Oscar.SubQuo(A,B)) # could go wrong
	end
end

@testset "Test simplify" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	A1 = R[x*y R(0)]
	B1 = R[R(0) R(1)]
	M1 = Oscar.SubQuo(A1,B1)
	M2,i2,p2 = Oscar.simplify(M1)
	#display(M1)
	#display(M2)
	for k=1:5
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:1])), M1)
		@test elem == i2(p2(elem))
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
		#elem = M2(R[randpoly(R) for i=1:ngens(M2)])
		@test elem == p2(i2(elem))
	end

	@test Oscar.iswelldefined(i2)
	@test Oscar.iswelldefined(p2)
	@test Oscar.isbijective(i2)
	@test Oscar.isbijective(p2)
	#@test inv(i2) === p2 && inv(p2) === i2

	A1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
	B1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
	M1 = Oscar.SubQuo(A1,B1)
	M2,i2,p2 = Oscar.simplify(M1)
	#display(M1)
	#display(M2)
	for k=1:5
		#elem = M1([randpoly(R) for i=1:3])
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:3])), M1)
		@test elem == i2(p2(elem))
		#elem = M2([randpoly(R) for i=1:ngens(M2)])
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test Oscar.iswelldefined(i2)
	@test Oscar.iswelldefined(p2)
	@test Oscar.isbijective(i2)
	@test Oscar.isbijective(p2)
	#@test inv(i2) === p2 && inv(p2) === i2

	A1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:3])
	B1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:3])
	M1 = Oscar.SubQuo(A1,B1)
	M2,i2,p2 = Oscar.simplify(M1)
	#display(M1)
	#display(M2)
	#println(matrix(p2))

	for k=1:5
		#elem = M1(R[randpoly(R) for i=1:3])
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1, i=1:3])), M1)
		@test elem == i2(p2(elem))
		#elem = M2(R[randpoly(R) for i=1:ngens(M2)])
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test Oscar.iswelldefined(i2)
	@test Oscar.iswelldefined(p2)
	@test Oscar.isbijective(i2)
	@test Oscar.isbijective(p2)
	#@test inv(i2) === p2 && inv(p2) === i2

	M1 = Oscar.SubQuo(B1,A1)
	M2,i2,p2 = Oscar.simplify(M1)
	#display(M1)
	#display(M2)
	for k=1:5
		#elem = M1(R[randpoly(R) for i=1:2])
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1, i=1:2])), M1)
		@test elem == i2(p2(elem))
		#elem = M2(R[randpoly(R) for i=1:ngens(M2)])
		elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test Oscar.iswelldefined(i2)
	@test Oscar.iswelldefined(p2)
	@test Oscar.isbijective(i2)
	@test Oscar.isbijective(p2)
	#@test inv(i2) === p2 && inv(p2) === i2
end

@testset "testing quotient modules" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])

  M1 = Oscar.SubQuo(R[x^2*y^3-x*y y^3 x^2*y; 2*x^2 3*y^2*x 4],R[x^4*y^5 x*y y^4])
  N1 = Oscar.SubQuo(R[x^4*y^5-4*x^2 x*y-6*y^2*x y^4-8],R[x^4*y^5 x*y y^4])
  Q1,p1 = quo(M1,N1,:store)

  @test Q1 == Oscar.SubQuo(R[x^2*y^3-x*y y^3 x^2*y],R[x^4*y^5 x*y y^4; x^4*y^5-4*x^2  -6*x*y^2+x*y  y^4-8])
  for k=1:5
    elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:2])), M1)  
    @test p1(elem) == Oscar.map_canonically(Q1, elem)
    #@test p1(elem) == Q1(elem)
  end

  M2 = Oscar.SubQuo(R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],R[x^3-y^2 y^4-x-y])
  elems = [Oscar.SubQuoElem(Oscar.sparse_row(R[x*y -x*y^2 x*y]), M2), Oscar.SubQuoElem(Oscar.sparse_row(R[x R(0) R(-1)]), M2)]
  Q2,p2 = quo(M2,elems,:store)

  @test Q2 == Oscar.SubQuo(R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],
            R[x^3-y^2 y^4-x-y; x^2*y^3+x^2*y^2+x^3*y^3+x*y^3-x^5*y^2 x^4*y+2*x*y^2-x*y^5+x^2*y^2; x^2*y-y^2 x^4+x*y])
  for k=1:5
    elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:3])), M2)
    @test p2(elem) == Oscar.map_canonically(Q2, elem)
  end

  M3 = Oscar.SubQuo(R[x^2*y+13*x*y+2x-1 x^4 2*x*y; y^4 3*x -1],R[y^2 x^3 y^2])
  N3 = Oscar.SubQuo(R[x^2*y+13*x*y+2x-1-x*y^2 0 x^4-x*y^2; y^4-x*y^2 3*x-x^4 -1-x*y^2],R[2*y^2 2*x^3 2*y^2])
  N3 = Oscar.SubQuo(R[x^2*y+13*x*y+2x-1-x*y^2 0 2*x*y-x*y^2; y^4-x*y^2 3*x-x^4 -1-x*y^2],R[2*y^2 2*x^3 2*y^2])
  Q3,p3 = quo(M3,N3,:store)

  @test iszero(quo(M3,M3))
  @test iszero(Q3)
  for k=1:5
    elem = Oscar.SubQuoElem(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:1])), M3)
    @test p3(elem) == Oscar.map_canonically(Q3, elem)
    @test iszero(p3(elem))
  end
end

@testset "testing submodules" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])

  M1 = Oscar.SubQuo(R[x^2*y+x*y x*y^2-x; x+x*y^2 y^3],R[x^2 y^3-x])
  S1,i1 = Oscar.sub(M1, [M1(Oscar.sparse_row(R[1 1])),M1(Oscar.sparse_row(R[y -x]))], :store)

  @test S1 == Oscar.SubQuo(R[x*y^2+x^3-x^2 x*y^3-x*y-x^2; x^2*y+x*y^2+x*y-x^2+x x*y^2],R[x^2 y^3-x])
  for k=1:5
      elem = S1(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:2])))
      @test i1(elem) == Oscar.map_canonically(M1, elem)
  end

  M2 = Oscar.SubQuo(R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],R[x^3-y^2 y^4-x-y])
  S2,i2 = sub(M2,[M2(Oscar.sparse_row(R[x*y -x*y^2 x*y])),M2(Oscar.sparse_row(R[x 0 -1]))], :store)

  @test S2 == Oscar.SubQuo(R[x^2*y^3+x^2*y^2+x^3*y^3+x*y^3-x^5*y^2 x^4*y+2*x*y^2-x*y^5+x^2*y^2; x^2*y-y^2 x^4+x*y],R[x^3-y^2 y^4-x-y])
  for k=1:5
      elem = S2(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1,i=1:2])))
      @test i2(elem) == Oscar.map_canonically(M2, elem)
  end

  M3 = Oscar.SubQuo(R[x*y^2 x^3+2*y; x^4 y^3; x*y+y^2 x*y],R[x^3-y^2 y^4-x-y])
  elems = [M3(Oscar.sparse_row(R[0 6 0])),M3(Oscar.sparse_row(R[9 0 -x])),M3(Oscar.sparse_row(R[0 0 -42]))]
  S3,i3 = sub(M3,elems,:store)

  @test S3 == M3
  for k=1:5
      elem = S3(Oscar.sparse_row(array_to_matrix([randpoly(R) for _=1:1, i=1:3])))
      @test i3(elem) == Oscar.map_canonically(M3, elem)
  end
end

@testset "testing Hom" begin
	R, (x0,x1,x2,x3,x4,x5) = PolynomialRing(QQ, ["x0", "x1", "x2", "x3", "x4", "x5"])
	f1=R[-x2*x3 -x4*x5 0; x0*x1 0 -x4*x5; 0 x0*x1 -x2*x3]'
	g1=R[x0*x1 x2*x3 x4*x5]'
	M = cokernel(f1)
	N = cokernel(g1)
	SQ = Oscar.hom(M,N)[1]

	#println("-------")
	#display(SQ)
	#println("--------")
	#display(simplify(SQ)[1])
	#println("--------")

	f2 = R[0 0]
	M = cokernel(f1)
	N = cokernel(f2)
	SQ = hom(M,N)[1]

	#println("-------")
	#display(SQ)
	#println("--------")
	#display(simplify(SQ)[1])
	#println("--------")

	function free_module_SQ(ring,n)
		A = one(AbstractAlgebra.MatrixSpace(ring,n,n))
		B = zero(AbstractAlgebra.MatrixSpace(ring,1,n))
		return Oscar.SubQuo(A,B)
	end

	# test if Hom(R^n,R^m) gives R^(m*n): (there is a dedicated function for free modules but this is a simple test for the function for subquos)
	for n=1:5
		for m=1:5
			M=free_module_SQ(R,n)
			N=free_module_SQ(R,m)
			#SQ = Oscar.hom(M,N)[1] # it's weird that it fails without simplifying
			SQ = Oscar.simplify(Oscar.hom(M,N)[1])[1]
			@test SQ==free_module_SQ(R,n*m)
		end
	end
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	A1 = R[x^2+1 x*y; x^2+y^3 x*y]
	B1 = R[x+x^4+y^2+1 x^5; y^4-3 x*y^2-1]
	M1 = Oscar.SubQuo(A1, B1)
	A2 = R[x;]
	B2 = R[y^3;]
	M2 = Oscar.SubQuo(A2,B2)
	SQ = Oscar.hom(M1,M2)[1]
	for v in Oscar.gens(SQ)
		@test v == Oscar.homomorphism_to_module_elem(SQ, Oscar.homomorphism(v))
	end

	# test if hom(zero-module, ...) is zero
	Z = Oscar.FreeMod(R,0)
	@test iszero(Oscar.hom(Z,Z)[1])
	for k=1:10
		A = array_to_matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
		B = array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		N = Oscar.SubQuo(A,B)

		@test iszero(hom(N,Z)[1])
		@test iszero(hom(Z,N)[1])
	end

	# test welldefinedness of randomly generated homomorphisms (using hom() and homomorphism())
	for k=1:10
		A1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
		A2 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:2])
		B1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		B2 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		N = Oscar.SubQuo(A1,B1)
		M = Oscar.SubQuo(A2,B2)
		HomNM = Oscar.hom(N,M)[1]
		for l=1:10
			H = HomNM(Oscar.sparse_row(array_to_matrix([randpoly(R,0:15,2,1) for _=1:1, j=1:AbstractAlgebra.ngens(HomNM)])))
			H = Oscar.homomorphism(H)
			@test Oscar.iswelldefined(H)
		end
	end
end

# testing lift ?

@testset "homomorphism testing" begin
	# This test doesn't terminate due to a bug in syz
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	k_max = 10
	for k=1:k_max
		print("\rhomomorphism testing: test ",k," out of ",k_max,"     ")
		A1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
		A2 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:2])
		B1 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		B2 = array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])


		#1) H: N --> M where N is a cokernel, H should be an isomorphism
		M = Oscar.SubQuo(A1,B1)
		N, H = Oscar.present_as_cokernel(M, :store)
		Hinv = H.inverse_isomorphism
		@test Oscar.iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = Oscar.kernel(H)
		ImH,iImH = Oscar.image(H)

		NmodKerH, pNmodKerH = quo(N,KerH, :store)
		Hbar = Oscar.SubQuoHom(NmodKerH,M,matrix(H))
		#Hbar = AbstractAlgebra.Generic.ModuleHomomorphism(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = Oscar.restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test Oscar.iswelldefined(Hbar)
		@test Oscar.isbijective(Hbar)

		Hbar_inv = Oscar.inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in Oscar.gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in Oscar.gens(ImH)])
		#######################################################################

		# test, if H is bijective with inverse Hinv:
		@test Oscar.isbijective(H)
		@test all([Oscar.inv(H)(H(g))==g for g in gens(N)])
		@test all([H(Oscar.inv(H)(g))==g for g in gens(M)])
		@test Oscar.inv(H) === Hinv
		@test ImH == M
		@test Oscar.iszero(KerH)


		#2) H: N --> M = N/(submodule of N) canonical projection
		M,H = Oscar.quo(N,[N(Oscar.sparse_row(R[1 x^2-1 x*y^2])),N(Oscar.sparse_row(R[y^3 y*x^2 x^3]))],:store)
		@test Oscar.iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = Oscar.kernel(H)
		ImH,iImH = Oscar.image(H)

		NmodKerH, pNmodKerH = Oscar.quo(N,KerH, :store)
		Hbar = Oscar.SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		#Hbar = AbstractAlgebra.Generic.ModuleHomomorphism(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = Oscar.restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test Oscar.iswelldefined(Hbar)
		@test Oscar.isbijective(Hbar)

		Hbar_inv = Oscar.inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
		#######################################################################

		# test, if H is surjective:
		@test ImH == M


		#3) H:N --> M neither injective nor surjective, also: tests 'restrict_domain()'
		MM = Oscar.SubQuo(A1,B1)
		M, iM, iSQ = Oscar.sum(MM, Oscar.SubQuo(A2,B1))
		NN, p, i = Oscar.direct_product(MM,Oscar.SubQuo(A2,B2), task = :both)
		i1,i2 = i[1],i[2]
		p1,p2 = p[1],p[2]
		nn = Oscar.ngens(NN)
		N,iN = Oscar.sub(NN,[NN(Oscar.sparse_row(array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:nn]))), NN(Oscar.sparse_row(array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:nn]))), NN(Oscar.sparse_row(array_to_matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:nn])))], :store)
		
		H = Oscar.restrict_domain(p1*iM,N)
		@test Oscar.iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = Oscar.kernel(H)
		ImH,iImH = Oscar.image(H)

		NmodKerH, pNmodKerH = Oscar.quo(N,KerH, :store)
		Hbar = Oscar.SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = Oscar.restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test Oscar.iswelldefined(Hbar)
		@test Oscar.isbijective(Hbar)

		Hbar_inv = Oscar.inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in Oscar.gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in Oscar.gens(ImH)])
		#######################################################################



		#4) H: M --> N random map created via the hom() function
		N = Oscar.SubQuo(A1,B1)
		M = Oscar.SubQuo(A2,B2)
		HomNM = Oscar.hom(N,M)[1]
		H = HomNM(Oscar.sparse_row(array_to_matrix([randpoly(R,0:15,2,1) for _=1:1,j=1:Oscar.ngens(HomNM)])))
		H = Oscar.homomorphism(H)
		@test Oscar.iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = Oscar.kernel(H)
		ImH,iImH = Oscar.image(H)

		NmodKerH, pNmodKerH = Oscar.quo(N,KerH, :store)
		Hbar = Oscar.SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		#Hbar = AbstractAlgebra.Generic.ModuleHomomorphism(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = Oscar.restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test Oscar.iswelldefined(Hbar)
		@test Oscar.isbijective(Hbar)

		Hbar_inv = Oscar.inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in Oscar.gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in Oscar.gens(ImH)])
		#######################################################################
	end
	print("\r                                        ")
	print("\r")
end