using Random

RNG = Random.MersenneTwister(42)

"""
	randpoly(R::Ring,coeffs=0:9,max_exp=4,max_terms=8)
> Return a random Polynomial from the Polynomial Ring `R` with coefficients in `coeffs`
> with exponents between `0` and `max_exp` und between `0` and `max_terms` terms
"""
function randpoly(R::Oscar.Ring,coeffs=0:9,max_exp=4,max_terms=8)
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
	matrix(A::Array,R::Ring = parent(A[1,1]))
Return `A` as an AbstractAlgebra Matrix
"""
function Oscar.matrix(A::Array,R::AbstractAlgebra.Ring = parent(A[1,1]))
	Mat = AbstractAlgebra.MatrixSpace(R,size(A)...)
	return Mat(R.(A))
end

@testset "Intersection of modules" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

  A1 = R[x y;
        2*x^2 3*y^2]
  A2 = R[x^3 x^2*y;
        (2*x^2+x) (2*y^3+y)]
  B = R[4*x*y^3 (2*x+y)]
  F2 = FreeMod(R,2)

  M1 = SubQuo(F2, A1, B)
  M2 = SubQuo(F2, A2, B)

  A1 = R[x R(1); y x^2]
  A2 = R[x y]
  B1 = R[x^2 y^2]
  res = R[(x*y^2) (y^3-x*y+y^2); (-x*y) (-x^2*y^2+x*y^3-x^2*y+y^3-x*y)]

  @test M1 == intersect(M1, M1)[1]
  @test SubQuo(F2, res, B1) == intersect(SubQuo(F2, A1, B1), SubQuo(F2, A2, B1))[1] #macaulay2

  A1 = R[x y; x^2 y^2]
  A2 = R[x+x^2 y+y^2]
  @test SubQuo(F2, A2, B1) == intersect(SubQuo(F2, A1,B1), SubQuo(F2, A2,B1))[1]
end

@testset "Presentation" begin

	# over Integers
	R, (x,y,z) = PolynomialRing(ZZ, ["x", "y", "z"])
	generator_matrices = [R[x x^2*y; y y^2*x^2], R[x x^2*y; y y^2*x^2], R[x y; x^2*y y^2*x], R[x+R(1) x^10; x^2+y^2 y^4-x^4], R[42*x*y 7*x^2; 6*x 9*y^2]]
	relation_matrices  = [R[x^3 y^4], R[x^3 y^4; x^2*y x*y^2], R[x*y^2 x; y^3 x*y^3], R[x x*y], R[3*x*y 7*x*y; 42 7]]
	true_pres_matrices = [R[x^5*y-y^4 -x^5+x*y^3], R[-x^2*y-x*y-y x^2+x; x*y^4-x^3*y+y^4-x^2*y -x*y^3+x^3; -y^4+x^2*y x*y^3-x^3; x^2*y^3+x*y^3+y^3 -x^2*y^2-x*y^2], R[-x*y R(1); -x^2*y^5+x*y^3 R(0)], R[x^5-x*y^4+x^3*y+x*y^3 x^11-x^2*y-x*y; -x^5*y^7+x*y^11+2*x^5*y^6-x^3*y^8-3*x*y^10+2*x^5*y^5+2*x^3*y^7-3*x^5*y^4+2*x^3*y^6+5*x*y^8+2*x^5*y^3-3*x^3*y^5-5*x*y^7+2*x^3*y^4+2*x*y^6 -x^11*y^7+2*x^11*y^6+2*x^11*y^5-3*x^11*y^4+2*x^11*y^3+x^2*y^8-2*x^2*y^7+x*y^8-2*x^2*y^6-2*x*y^7+3*x^2*y^5-2*x*y^6-2*x^2*y^4+3*x*y^5-2*x*y^4], R[-13*y R(0); -377 -2639*x; -39 -273*x; -13*y-39 -273*x; -13 -91*x; y^2-42*x -294*x^2+21*x*y; 9*y^2-x+26*y+78 -7*x^2+189*x*y+546*x; -y^2+3*x 21*x^2-21*x*y]]
	for (A,B,true_pres_mat) in zip(generator_matrices, relation_matrices, true_pres_matrices)
		SQ = SubQuo(A,B)
		pres_mat = matrix(present_as_cokernel(SQ).quo)
		F = FreeMod(R,ncols(pres_mat))
		@test cokernel(F,pres_mat) == cokernel(F,true_pres_mat)

		pres_SQ, i = present_as_cokernel(SQ, :both)
		p = i.inverse_isomorphism
		@test iswelldefined(i)
		@test iswelldefined(p)
		@test isbijective(i)
		@test isbijective(p)
	end


	# over Rationals
	R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	generator_matrices = [R[x x^2*y; y y^2*x^2], R[x x^2*y; y y^2*x^2], R[x y; x^2*y y^2*x], R[x+R(1) x^10; x^2+y^2 y^4-x^4], R[x+R(1) x^10; x^2+y^2 y^4-x^4]]
	relation_matrices  = [R[x^3 y^4], R[x^3 y^4; x^2*y x*y^2], R[x*y^2 x; y^3 x*y^3], R[x+y x+y; x*y x*y], R[x+y x+y; y x; x y]]
	true_pres_matrices = [R[x^5*y-y^4 -x^5+x*y^3], R[-x^2*y-x*y-y x^2+x; -x^2*y^4+x^4*y x^2*y^3-x^4], R[-x*y R(1); -x^2*y^5+x*y^3 R(0)], R[-x^4+y^4-x^2-y^2 -x^10+x+R(1)], R[R(0) -x^2+y^2; -2*x^2 -x^9*y+x+1; -2*x*y^2 -x^8*y^3+y^2+x; -2*x*y^2+2*y^2 -x^8*y^3+x^7*y^3+y^2-1]]
	for (A,B,true_pres_mat) in zip(generator_matrices, relation_matrices, true_pres_matrices)
		SQ = SubQuo(A,B)
		pres_mat = matrix(present_as_cokernel(SQ).quo)
		F = FreeMod(R,ncols(pres_mat))
		@test cokernel(F,pres_mat) == cokernel(F,true_pres_mat)

		pres_SQ, i = present_as_cokernel(SQ, :both)
    	p = i.inverse_isomorphism
		@test iswelldefined(i)
		@test iswelldefined(p)
		@test isbijective(i)
		@test isbijective(p)
	end
end


@testset "Test kernel" begin

	# over Integers
	R, (x,y,z) = PolynomialRing(ZZ, ["x", "y", "z"])
	matrices = [R[x^2+x y^2+y; x^2+y y^2; x y], R[5*x^5+x*y^2 4*x*y+y^2+R(1); 4*x^2*y 2*x^2+3*y^2-R(5)]]
	kernels = [R[x^2-x*y+y -x^2+x*y -x^2+x*y-y^2-y], R[R(0) R(0)]]
	for (A,Ker) in zip(matrices, kernels)
		F1 = FreeMod(R, nrows(A))
		F2 = FreeMod(R, ncols(A))
		K,emb = kernel(FreeModuleHom(F1,F2,A))
    	@test K == image(emb)[1]
    	@test image(emb)[1] == image(map(F1,Ker))[1]
	end
	for k=1:3
		A = matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
    	A = map(A)
		K,emb = kernel(A)
		@test iszero(emb*A)
	end

	# over Rationals
	R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	matrices = [R[x^2+x y^2+y; x^2+y y^2; x y], R[5*x^5+x*y^2 4*x*y+y^2+R(1); 4*x^2*y 2*x^2+3*y^2-R(5)],
				R[8*x^2*y^2*z^2+13*x*y*z^2  12*x^2+7*y^2*z;
				13*x*y^2+12*y*z^2  4*x^2*y^2*z+8*x*y*z;
				9*x*y^2+4*z  12*x^2*y*z^2+9*x*y^2*z]]
	kernels = [R[x^2-x*y+y -x^2+x*y -x^2+x*y-y^2-y], R[R(0) R(0)], 
			   R[-36*x^3*y^4*z+156*x^3*y^3*z^2+144*x^2*y^2*z^4+117*x^2*y^4*z+108*x*y^3*z^3-72*x^2*y^3*z-16*x^2*y^2*z^2-32*x*y*z^2 -96*x^4*y^3*z^4-72*x^3*y^4*z^3-156*x^3*y^2*z^4-117*x^2*y^3*z^3+63*x*y^4*z+108*x^3*y^2+28*y^2*z^2+48*x^2*z 32*x^4*y^4*z^3+116*x^3*y^3*z^3+104*x^2*y^2*z^3-91*x*y^4*z-84*y^3*z^3-156*x^3*y^2-144*x^2*y*z^2]]
	for (A,Ker) in zip(matrices, kernels)
		F1 = FreeMod(R, nrows(A))
		F2 = FreeMod(R, ncols(A))
		K,emb = kernel(FreeModuleHom(F1,F2,A))
    	@test K == image(emb)[1]
    	@test image(emb)[1] == image(map(F1,Ker))[1]
	end
	for k=1:3
		A = matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
    	A = map(A)
		K,emb = kernel(A)
    	@test image(emb)[1] == K
		@test iszero(emb*A)
	end
end

@testset "iszero(SubQuo)" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	A = R[x^2+2*x*y y^2*x-2*x^2*y;-y x*y]
	B = R[x^2 y^2*x;-y x*y]
	@test iszero(SubQuo(A,B))
	for k=1:3
		A = matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
		B = matrix([randpoly(R,0:15,2,2) for i=1:2,j=1:2])
		@test iszero(SubQuo(A,A))
		@test !iszero(SubQuo(A,B)) # could go wrong
	end
end

@testset "simplify subquotient" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	A1 = R[x*y R(0)]
	B1 = R[R(0) R(1)]
	M1 = SubQuo(A1,B1)
	M2,i2,p2 = simplify(M1)
	for k=1:5
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:1])), M1)
		@test elem == i2(p2(elem))
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test iswelldefined(i2)
	@test iswelldefined(p2)
	@test isbijective(i2)
	@test isbijective(p2)

	A1 = matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
	B1 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
	M1 = SubQuo(A1,B1)
	M2,i2,p2 = simplify(M1)
	for k=1:5
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:3])), M1)
		@test elem == i2(p2(elem))
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test iswelldefined(i2)
	@test iswelldefined(p2)
	@test isbijective(i2)
	@test isbijective(p2)

	A1 = matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:3])
	B1 = matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:3])
	M1 = SubQuo(A1,B1)
	M2,i2,p2 = simplify(M1)

	for k=1:5
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:3])), M1)
		@test elem == i2(p2(elem))
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test iswelldefined(i2)
	@test iswelldefined(p2)
	@test isbijective(i2)
	@test isbijective(p2)

	M1 = SubQuo(B1,A1)
	M2,i2,p2 = simplify(M1)
	for k=1:5
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:2])), M1)
		@test elem == i2(p2(elem))
		elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
		@test elem == p2(i2(elem))
	end

	@test iswelldefined(i2)
	@test iswelldefined(p2)
	@test isbijective(i2)
	@test isbijective(p2)
end

@testset "quotient modules" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])

  F3 = FreeMod(R,3)
  M1 = SubQuo(F3,R[x^2*y^3-x*y y^3 x^2*y; 2*x^2 3*y^2*x 4],R[x^4*y^5 x*y y^4])
  N1 = SubQuo(F3,R[x^4*y^5-4*x^2 x*y-6*y^2*x y^4-8],R[x^4*y^5 x*y y^4])
  Q1,p1 = quo(M1,N1,:store)

  @test Q1 == SubQuo(F3,R[x^2*y^3-x*y y^3 x^2*y],R[x^4*y^5 x*y y^4; x^4*y^5-4*x^2  -6*x*y^2+x*y  y^4-8])
  for k=1:5
    elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:2])), M1)  
    @test p1(elem) == map_canonically(Q1, elem)
  end

  F2 = FreeMod(R,2)
  M2 = SubQuo(F2,R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],R[x^3-y^2 y^4-x-y])
  elems = [SubQuoElem(sparse_row(R[x*y -x*y^2 x*y]), M2), SubQuoElem(sparse_row(R[x R(0) R(-1)]), M2)]
  Q2,p2 = quo(M2,elems,:store)

  @test Q2 == SubQuo(F2,R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],
            R[x^3-y^2 y^4-x-y; x^2*y^3+x^2*y^2+x^3*y^3+x*y^3-x^5*y^2 x^4*y+2*x*y^2-x*y^5+x^2*y^2; x^2*y-y^2 x^4+x*y])
  for k=1:5
    elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:3])), M2)
    @test p2(elem) == map_canonically(Q2, elem)
  end

  M3 = SubQuo(F3,R[x^2*y+13*x*y+2x-1 x^4 2*x*y; y^4 3*x -1],R[y^2 x^3 y^2])
  #N3 = SubQuo(F3,R[x^2*y+13*x*y+2x-1-x*y^2 0 x^4-x*y^2; y^4-x*y^2 3*x-x^4 -1-x*y^2],R[2*y^2 2*x^3 2*y^2])
  N3 = SubQuo(F3,R[x^2*y+13*x*y+2x-1-x*y^2 0 2*x*y-x*y^2; y^4-x*y^2 3*x-x^4 -1-x*y^2],R[2*y^2 2*x^3 2*y^2])
  Q3,p3 = quo(M3,N3,:store)

  @test iszero(quo(M3,M3))
  @test iszero(Q3)
  for k=1:5
    elem = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:1])), M3)
    @test p3(elem) == map_canonically(Q3, elem)
    @test iszero(p3(elem))
  end
end

@testset "submodules" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])

  F2 = FreeMod(R,2)
  M1 = SubQuo(F2,R[x^2*y+x*y x*y^2-x; x+x*y^2 y^3],R[x^2 y^3-x])
  S1,i1 = sub(M1, [M1(sparse_row(R[1 1])),M1(sparse_row(R[y -x]))], :store)

  @test S1 == SubQuo(F2,R[x*y^2+x^3-x^2 x*y^3-x*y-x^2; x^2*y+x*y^2+x*y-x^2+x x*y^2],R[x^2 y^3-x])
  for k=1:5
      elem = S1(sparse_row(matrix([randpoly(R) for _=1:1,i=1:2])))
      @test i1(elem) == map_canonically(M1, elem)
  end

  M2 = SubQuo(F2,R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],R[x^3-y^2 y^4-x-y])
  S2,i2 = sub(M2,[M2(sparse_row(R[x*y -x*y^2 x*y])),M2(sparse_row(R[x 0 -1]))], :store)

  @test S2 == SubQuo(F2,R[x^2*y^3+x^2*y^2+x^3*y^3+x*y^3-x^5*y^2 x^4*y+2*x*y^2-x*y^5+x^2*y^2; x^2*y-y^2 x^4+x*y],R[x^3-y^2 y^4-x-y])
  for k=1:5
      elem = S2(sparse_row(matrix([randpoly(R) for _=1:1,i=1:2])))
      @test i2(elem) == map_canonically(M2, elem)
  end

  M3 = SubQuo(F2,R[x*y^2 x^3+2*y; x^4 y^3; x*y+y^2 x*y],R[x^3-y^2 y^4-x-y])
  elems = [M3(sparse_row(R[0 6 0])),M3(sparse_row(R[9 0 -x])),M3(sparse_row(R[0 0 -42]))]
  S3,i3 = sub(M3,elems,:store)

  @test S3 == M3
  for k=1:5
      elem = S3(sparse_row(matrix([randpoly(R) for _=1:1, i=1:3])))
      @test i3(elem) == map_canonically(M3, elem)
  end
end

@testset "Hom module" begin
	R, (x0,x1,x2,x3,x4,x5) = PolynomialRing(QQ, ["x0", "x1", "x2", "x3", "x4", "x5"])
	f1=R[-x2*x3 -x4*x5 0; x0*x1 0 -x4*x5; 0 x0*x1 -x2*x3]'
	g1=R[x0*x1 x2*x3 x4*x5]'
	M = cokernel(f1)
	N = cokernel(g1)
	SQ = hom(M,N)[1]

	f2 = R[0 0]
	M = cokernel(f1)
	N = cokernel(f2)
	SQ = hom(M,N)[1]


	function free_module_SQ(ring,n)
		F = FreeMod(ring,n)
		return SubQuo(F,gens(F))
	end

	function free_module_SQ(F::FreeMod)
		return SubQuo(F,gens(F))
	end

	# test if Hom(R^n,R^m) gives R^(m*n): (there is a dedicated function for free modules but this is a simple test for the function for subquos)
	for n=1:5
		for m=1:5
			M=free_module_SQ(R,n)
			N=free_module_SQ(R,m)
			SQ = hom(M,N)[1]
			#SQ = simplify(hom(M,N)[1])[1]
			@test SQ == free_module_SQ(free_module(SQ))
		end
	end
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	A1 = R[x^2+1 x*y; x^2+y^3 x*y]
	B1 = R[x+x^4+y^2+1 x^5; y^4-3 x*y^2-1]
	M1 = SubQuo(A1, B1)
	A2 = R[x;]
	B2 = R[y^3;]
	M2 = SubQuo(A2,B2)
	SQ = hom(M1,M2)[1]
	for v in gens(SQ)
		@test v == module_elem(SQ, homomorphism(v))
	end

	# test if hom(zero-module, ...) is zero
	Z = FreeMod(R,0)
	@test iszero(hom(Z,Z)[1])
	for k=1:10
		A = matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
		B = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		N = SubQuo(A,B)

		@test iszero(hom(N,Z)[1])
		@test iszero(hom(Z,N)[1])
	end

	# test welldefinedness of randomly generated homomorphisms (using hom() and homomorphism())
	for k=1:10
		A1 = matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
		A2 = matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:2])
		B1 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		B2 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		N = SubQuo(A1,B1)
		M = SubQuo(A2,B2)
		HomNM = hom(N,M)[1]
		for l=1:10
			H = HomNM(sparse_row(matrix([randpoly(R,0:15,2,1) for _=1:1, j=1:AbstractAlgebra.ngens(HomNM)])))
			H = homomorphism(H)
			@test iswelldefined(H)
		end
	end
end

@testset "tensoring morphisms" begin
	R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])

	F2 = FreeMod(R,2)
	F3 = FreeMod(R,3)
	F4 = FreeMod(R,4)

	for _=1:10
		A1 = matrix([randpoly(R,0:15,4,3) for i=1:3,j=1:2])
		B1 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])

		
		A2 = matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:3])
		B2 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:3])

		M1 = SubQuo(F2,A1,B1)
		M2 = SubQuo(F3,A2,B2)

		M,pure_M = tensor_product(M1,M2, task=:map)
		phi = hom_tensor(M,M,[identity_map(M1),identity_map(M2)])

		for _=1:3
			v = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M)])), M)
			@test phi(v) == v
		end

		A3 = matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:2])
		#B3 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		M3 = SubQuo(Oscar.SubModuleOfFreeModule(F2,A3))

		N,pure_N = tensor_product(M3,F4, task=:map)

		M3_to_M1 = SubQuoHom(M3,M1, matrix([randpoly(R,0:2,2,2) for i=1:ngens(M3), j=1:ngens(M1)]))
		@assert iswelldefined(M3_to_M1)
		F4_to_M2 = FreeModuleHom(F4,M2, matrix([randpoly(R,0:2,2,2) for i=1:ngens(F4), j=1:ngens(M2)]))

		phi = hom_tensor(N,M,[M3_to_M1,F4_to_M2])
		u1 = SubQuoElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M3)])), M3)
		u2 = FreeModElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(F4)])), F4)
		@test phi(pure_N((u1,u2))) == pure_M((M3_to_M1(u1),F4_to_M2(u2)))
	end
end

@testset "direct product" begin
	R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	
	F2 = FreeMod(R,2)
	F3 = FreeMod(R,3)

	for _=1:3
		A1 = matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
		B1 = matrix([randpoly(R,0:15,2,2) for i=1:1,j=1:2])
		M1 = SubQuo(F2,A1,B1)

		A2 = matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:3])
		B2 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:3])
		M2 = SubQuo(F3,A2,B2)

		prod_M, proj, emb = direct_product(M1,M2,task=:both)
		@test length(proj) == length(emb) == 2
		@test ngens(prod_M) == ngens(M1) + ngens(M2)

		for g in gens(prod_M)
			@test g == sum([emb[i](proj[i](g)) for i=1:length(proj)])
		end
		for g in gens(M1)
			@test g == proj[1](emb[1](g))
		end
		for g in gens(M2)
			@test g == proj[2](emb[2](g))
		end

		A1 = matrix([randpoly(R,0:15,2,2) for i=1:3,j=1:2])
		B1 = matrix([randpoly(R,0:15,2,2) for i=1:1,j=1:2])
		N1 = SubQuo(F2,A1,B1)

		A2 = matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:3])
		B2 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:3])
		N2 = SubQuo(F3,A2,B2)

		prod_N = direct_product(N1,N2)
		@test ngens(prod_M) == ngens(M1) + ngens(M2)

		for g in gens(prod_N)
			@test g == sum([Hecke.canonical_injection(prod_N,i)(Hecke.canonical_projection(prod_N,i)(g)) for i=1:2])
		end
		for g in gens(N1)
			@test g == Hecke.canonical_projection(prod_N,1)(Hecke.canonical_injection(prod_N,1)(g))
		end
		for g in gens(N2)
			@test g == Hecke.canonical_projection(prod_N,2)(Hecke.canonical_injection(prod_N,2)(g))
		end

		# testing hom_prod_prod

		M1_to_N1 = SubQuoHom(M1,N1,zero_matrix(R,3,3))
		H12 = hom(M1,N2)[1]
		H21 = hom(M2,N1)[1]
		M1_to_N2 = iszero(H12) ? SubQuoHom(M1,N2,zero_matrix(R,3,2)) : homomorphism(H12[1])
		M2_to_N1 = iszero(H21) ? SubQuoHom(M2,N1,zero_matrix(R,2,3)) : homomorphism(H21[1])
		M2_to_N2 = SubQuoHom(M2,N2,R[0 0; 1 0])
		@assert iswelldefined(M1_to_N1)
		@assert iswelldefined(M1_to_N2)
		@assert iswelldefined(M2_to_N1)
		@assert iswelldefined(M2_to_N2)

		phi = hom_prod_prod(prod_M,prod_N,[M1_to_N1 M1_to_N2; M2_to_N1 M2_to_N2])
		for g in gens(M1)
			@test M1_to_N1(g) == Hecke.canonical_projection(prod_N,1)(phi(emb[1](g)))
			@test M1_to_N2(g) == Hecke.canonical_projection(prod_N,2)(phi(emb[1](g)))
		end
		for g in gens(M2)
			@test M2_to_N1(g) == Hecke.canonical_projection(prod_N,1)(phi(emb[2](g)))
			@test M2_to_N2(g) == Hecke.canonical_projection(prod_N,2)(phi(emb[2](g)))
		end
	end
end

@testset "Coordinates (lift)" begin
	Z3, a = FiniteField(3,1,"a")
	R, (x,y) = PolynomialRing(Z3, ["x", "y"])
	coeffs = [Z3(i) for i=0:1]

	A = R[x*y x^2+y^2; y^2 x*y;x^2+1 1]
	B = R[2*x^2 x+y; x^2+y^2-1 x^2+2*x*y]
	F = FreeMod(R,2)
	M = SubQuo(F,A,B)

	monomials = [x,y]
	coeff_generator = ([c1,c2] for c1 in coeffs for c2 in coeffs)
	for coefficients in ([c1,c2,c3] for c1 in coeff_generator for c2 in coeff_generator for c3 in coeff_generator)
		v = sparse_row(R, [(i,sum(coefficients[i][j]*monomials[j] for j=1:2)) for i=1:3])
		v_as_FreeModElem = sum([v[i]*repres(M[i]) for i=1:ngens(M)])

		elem1 = SubQuoElem(v_as_FreeModElem,M) 
		elem2 = SubQuoElem(v,M)

		@test elem1 == elem2
	end
end

@testset "module homomorphisms" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])

	F3 = FreeMod(R,3)
	F4 = FreeMod(R,4)
	phi = hom(F3,F4, [F4[1],F4[3]+x*F4[4],(x+y)*F4[4]] )
	@test iszero(preimage(phi,zero(F4)))
	@test phi(preimage(phi,y*(F4[3]+x*F4[4]+(x+y)*F4[3]))) == y*(F4[3]+x*F4[4]+(x+y)*F4[3])

	k_max = 8
	for k=1:k_max
		print("\rhomomorphism testing: test ",k," out of ",k_max,"     ")
		A1 = matrix([randpoly(R,0:15,2,1) for i=1:3,j=1:2])
		A2 = matrix([randpoly(R,0:15,2,1) for i=1:2,j=1:2])
		B1 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])
		B2 = matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:2])


		#1) H: N --> M where N is a cokernel, H should be an isomorphism
		M = SubQuo(A1,B1)
		N, H = present_as_cokernel(M, :store)
		Hinv = H.inverse_isomorphism
		@test iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = kernel(H)
		ImH,iImH = image(H)

		NmodKerH, pNmodKerH = quo(N,KerH, :store)
		Hbar = SubQuoHom(NmodKerH,M,matrix(H))
		Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test iswelldefined(Hbar)
		@test isbijective(Hbar)

		Hbar_inv = inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
		#######################################################################

		# test, if H is bijective with inverse Hinv:
		@test isbijective(H)
		@test all([inv(H)(H(g))==g for g in gens(N)])
		@test all([H(inv(H)(g))==g for g in gens(M)])
		@test inv(H) === Hinv
		@test ImH == M
		@test iszero(KerH)


		#2) H: N --> M = N/(submodule of N) canonical projection
		M,H = quo(N,[N(sparse_row(R[1 x^2-1 x*y^2])),N(sparse_row(R[y^3 y*x^2 x^3]))],:store)
		@test iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = kernel(H)
		ImH,iImH = image(H)

		NmodKerH, pNmodKerH = quo(N,KerH, :store)
		Hbar = SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test iswelldefined(Hbar)
		@test isbijective(Hbar)

		Hbar_inv = inv(Hbar)

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
		MM = SubQuo(A1,B1)
		M, iM, iSQ = sum(MM, SubQuo(A2,B1))
		NN, p, i = direct_product(MM,SubQuo(A2,B2), task = :both)
		i1,i2 = i[1],i[2]
		p1,p2 = p[1],p[2]
		nn = ngens(NN)
		N,iN = sub(NN,[NN(sparse_row(matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:nn]))), NN(sparse_row(matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:nn]))), NN(sparse_row(matrix([randpoly(R,0:15,2,1) for i=1:1,j=1:nn])))], :store)
		
		H = restrict_domain(p1*iM,N)
		@test iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = kernel(H)
		ImH,iImH = image(H)

		NmodKerH, pNmodKerH = quo(N,KerH, :store)
		Hbar = SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test iswelldefined(Hbar)
		@test isbijective(Hbar)

		Hbar_inv = inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
		#######################################################################



		#4) H: M --> N random map created via the hom() function
		N = SubQuo(A1,B1)
		M = SubQuo(A2,B2)
		HomNM = hom(N,M)[1]
		H = HomNM(sparse_row(matrix([randpoly(R,0:15,2,1) for _=1:1,j=1:ngens(HomNM)])))
		H = homomorphism(H)
		@test iswelldefined(H)

		## testing the homomorphism theorem: #################################
		KerH,iKerH = kernel(H)
		ImH,iImH = image(H)

		NmodKerH, pNmodKerH = quo(N,KerH, :store)
		Hbar = SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
		Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

		@test iswelldefined(Hbar)
		@test isbijective(Hbar)

		Hbar_inv = inv(Hbar)

		# test, if caching of inverse maps works:
		@test Hbar_inv === Hbar.inverse_isomorphism
		@test Hbar === Hbar_inv.inverse_isomorphism

		# test if Hbar and Hbar_inv are inverse to each other:
		@test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
		@test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
		#######################################################################
	end
	print("\r                                        ")
	print("\r")
end
