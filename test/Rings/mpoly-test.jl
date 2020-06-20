@testset "Polynomial Orderings" begin

	Qx, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	t = gen(Hecke.Globals.Qx)
	k1 , l= number_field(t^5+t^2+2)
	NFx = PolynomialRing(k1, ["x", "y", "z"])
	k2 = Nemo.GF(23)
	GFx = PolynomialRing(k2, ["x", "y", "z"])
	RNmodx=PolynomialRing(Nemo.ResidueRing(ZZ,17), :x => 1:2)[1]
	
	function rmIdeals(field :: Symbol)
		Ideals=[]
		if field == :Q
			Ox, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
			global Qx = Ox
			coeff = rand(-10:10,4,4)
			Oscar.push!(Ideals,ideal([3x^3*y+x^3+x*y^3+y^2*z^2,2x^3*z-x*y-x*z^3-y^4-z^2,2x^2*y*z-2x*y^2+x*z^2-y^4]))
		end
		if field == :NF
			Ox , (x, y, z) = PolynomialRing(k1, ["x", "y", "z"])
			coeff=fill(zero(k1),4,4)
			for k = 1:4
				for j= 1:4
					coeff[j,k]=dot(rand(-10:10,5),[one(k1),l,l^2,l^3,l^4])
				end
			end
			Oscar.push!(Ideals, ideal([2*x*y^4*z^2+(l-1)*x^2*y^3*z+(2*l)*x*y*z^2+7*y^3+(7*l+1),2*x^2*y^4*z+(l)*x^2*y*z^2-x*y^2*z^2+(2*l+3)*x^2*y*z-12*x+(12*l)*y,(2*l)*y^5*z+x^2*y^2*z-x*y^3*z+(-l)*x*y^3+y^4+2*y^2*z,(3*l)*x*y^4*z^3+(l+1)*x^2*y^2*z-x*y^3*z+4*y^3*z^2+(3*l)*x*y*z^3+4*z^2-x+(l)*y]))
		end
		if field == :GF
			Ox , (x,y,z) = PolynomialRing(k2, ["x", "y", "z"])
			init = rand(0:22,4,4)
			coeff=[(k2)(x) for x= init]
		end	
		for k= 1:4
			I=ideal([zero(Ox)])
			Oscar.groebner_assure(I)
			for j=1 : 4
				counter = 0
				e=rand(0:6,1,3)
				mon=coeff[j,k]*x^(e[1])*y^(e[2])*z^(e[3])
				red_mon = reduce(convert( I.gens.Sx, mon),I.gb.S)
				while red_mon == 0 && counter < 1000
					e=rand(0:6,1,3)
					mon=x^(e[1])*y^(e[2])*z^(e[3])
					red_mon = reduce(convert( I.gens.Sx, mon),I.gb.S)
					counter += 1
				end
				I=ideal([convert(Ox, x) for x= Oscar.push!(collect(gens(I.gb.S)), red_mon)])
				Oscar.groebner_assure(I)		
			end
			Oscar.push!(Ideals,I)
		end
		return Ideals
	end


	IdealsQx = rmIdeals(:Q)
	PolyQx = collect(gens(IdealsQx[2]))
	Poly = collect(gens(IdealsQx[1]))

	IdealsNFx = rmIdeals(:NF)
	PolyNFx = collect(gens(IdealsNFx[2]))

	IdealsGFx = rmIdeals(:GF)
	PolyGFx = collect(gens(IdealsGFx[2]))

	RFx,(x,y,z)= PolynomialRing(QQ,:x => 1:3)
	f=x^2*y+x*y^2+x*y+y^4

	@test Oscar._isless_lex(f,2,1)
	@test !Oscar._isless_lex(f,1,2)
	@test Oscar._isless_neglex(f,1,2)
	@test !Oscar._isless_neglex(f,2,1)
	@test Oscar._isless_revlex(f,3,1)
	@test !Oscar._isless_revlex(f,1,3)
	@test Oscar._isless_negrevlex(f,1,3)
	@test !Oscar._isless_negrevlex(f,3,1)
	@test Oscar._isless_deglex(f,1,4)
	@test !Oscar._isless_deglex(f,4,1)
	@test Oscar._isless_deglex(f,2,1)
	@test Oscar._isless_degrevlex(f,1,4)
	@test !Oscar._isless_degrevlex(f,4,1)
	@test Oscar._isless_degrevlex(f,2,1)
	@test !Oscar._isless_negdeglex(f,1,4)
	@test Oscar._isless_negdeglex(f,4,1)
	@test Oscar._isless_negdeglex(f,2,1)
	@test !Oscar._isless_negdeglex(f,1,4)
	@test Oscar._isless_negdegrevlex(f,4,1)
	@test Oscar._isless_negdegrevlex(f,2,1)
	@test Oscar.weighted_degree(f,2,[3,2,2]) == 7
	@test !Oscar._isless_weightlex(f,2,3,[3,2,2])
	@test Oscar._isless_weightlex(f,3,2,[3,2,2])
	@test !Oscar._isless_weightrevlex(f,2,3,[3,2,2])
	@test Oscar._isless_weightrevlex(f,3,2,[3,2,2])
	@test !Oscar._isless_weightneglex(f,2,3,[3,2,2])
	@test Oscar._isless_weightneglex(f,3,2,[3,2,2])
	@test Oscar._isless_weightnegrevlex(f,2,3,[3,2,2])
	@test !Oscar._isless_weightnegrevlex(f,3,2,[3,2,2])
	
	
	M=[0 0 0 ;1 1 1 ;1 2 3 ]
	@test !Oscar._isless_matrix(f,2,3,M)
	@test Oscar._isless_matrix(f,3,2,M)
	Oscar._perm_of_terms(f,Oscar._isless_deglex)

	Oscar.lt_from_ordering(RFx, :Wp, [3,2,2])
	Oscar.lt_from_ordering(RFx, :wp, [3,2,2])
	Oscar.lt_from_ordering(RFx, :Ws, [3,2,2])
	Oscar.lt_from_ordering(RFx, :ws, [3,2,2])

	(terms)(f, :lex)
	(coeffs)(f, :deglex)
	(exponent_vectors)(f,M)
	(monomials)(f, :weightlex, [3,2,2])


	Oscar.lt_from_ordering(RFx, :lex)
	Oscar.lt_from_ordering(RFx, :revlex)
	Oscar.lt_from_ordering(RFx, :deglex)
	Oscar.lt_from_ordering(RFx, :degrevlex)
	Oscar.lt_from_ordering(RFx, :neglex)
	Oscar.lt_from_ordering(RFx, :negrevlex)
	Oscar.lt_from_ordering(RFx, :negdeglex)
	Oscar.lt_from_ordering(RFx, :negdegrevlex)
	Oscar.lt_from_ordering(RFx,M)

	Oscar.terms(f,Oscar._isless_degrevlex)
	Oscar.coeffs(f,Oscar._isless_negdeglex)
	Oscar.exponent_vectors(f,Oscar._isless_revlex)
	Oscar.monomials(f,Oscar._isless_lex)
	
	A=[]
	Oscar.push!(A, PolyQx[1])

	@test Oscar.leading_term(f) == x^2*y
	@test Oscar.leading_coeff(f) == 1
	@test Oscar.leading_monomial(f) == x^2*y
	Oscar.leading_ideal([x for x= Poly])
	##Oscar.leading_ideal(PolyQx)
	Oscar.leading_ideal(A)
	Oscar.leading_ideal(parent(PolyQx[1]), A)
	Oscar.leading_ideal(IdealsQx[2])
	Oscar.leading_ideal(IdealsQx[2],:lex)
	factor(f)

end


@testset "Ideals" begin

	Qx, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	t = gen(Hecke.Globals.Qx)
	k1 , l= number_field(t^5+t^2+2)
	NFx = PolynomialRing(k1, ["x", "y", "z"])
	k2 = Nemo.GF(23)
	GFx = PolynomialRing(k2, ["x", "y", "z"])
	RNmodx=PolynomialRing(Nemo.ResidueRing(ZZ,17), :x => 1:2)[1]

	function rmIdeals(field :: Symbol)
		Ideals=[]
		if field == :Q
			Ox, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
			global Qx = Ox
			coeff = rand(-10:10,4,4)
			Oscar.push!(Ideals,ideal([3x^3*y+x^3+x*y^3+y^2*z^2,2x^3*z-x*y-x*z^3-y^4-z^2,2x^2*y*z-2x*y^2+x*z^2-y^4]))
		end
		if field == :NF
			Ox , (x, y, z) = PolynomialRing(k1, ["x", "y", "z"])
			coeff=fill(zero(k1),4,4)
			for k = 1:4
				for j= 1:4
					coeff[j,k]=dot(rand(-10:10,5),[one(k1),l,l^2,l^3,l^4])
				end
			end
			Oscar.push!(Ideals, ideal([2*x*y^4*z^2+(l-1)*x^2*y^3*z+(2*l)*x*y*z^2+7*y^3+(7*l+1),2*x^2*y^4*z+(l)*x^2*y*z^2-x*y^2*z^2+(2*l+3)*x^2*y*z-12*x+(12*l)*y,(2*l)*y^5*z+x^2*y^2*z-x*y^3*z+(-l)*x*y^3+y^4+2*y^2*z,(3*l)*x*y^4*z^3+(l+1)*x^2*y^2*z-x*y^3*z+4*y^3*z^2+(3*l)*x*y*z^3+4*z^2-x+(l)*y]))
		end
		if field == :GF
			Ox , (x,y,z) = PolynomialRing(k2, ["x", "y", "z"])
			init = rand(0:22,4,4)
			coeff=[(k2)(x) for x= init]
		end	
		for k= 1:4
			I=ideal([zero(Ox)])
			Oscar.groebner_assure(I)
			for j=1 : 4
				counter = 0
				e=rand(0:6,1,3)
				mon=coeff[j,k]*x^(e[1])*y^(e[2])*z^(e[3])
				red_mon = reduce(convert( I.gens.Sx, mon),I.gb.S)
				while red_mon == 0 && counter < 1000
					e=rand(0:6,1,3)
					mon=x^(e[1])*y^(e[2])*z^(e[3])
					red_mon = reduce(convert( I.gens.Sx, mon),I.gb.S)
					counter += 1
				end
				I=ideal([convert(Ox, x) for x= Oscar.push!(collect(gens(I.gb.S)), red_mon)])
				Oscar.groebner_assure(I)		
			end
			Oscar.push!(Ideals,I)
		end
		return Ideals
	end

	IdealsQx = rmIdeals(:Q)
	PolyQx = collect(gens(IdealsQx[2]))
	Poly = collect(gens(IdealsQx[1]))

	IdealsNFx = rmIdeals(:NF)
	PolyNFx = collect(gens(IdealsNFx[2]))

	IdealsGFx = rmIdeals(:GF)
	PolyGFx = collect(gens(IdealsGFx[2]))

	@test issubset(intersect(IdealsQx[3],IdealsQx[4] + IdealsQx[3]),(IdealsQx[1]+IdealsQx[2])*IdealsQx[3])
	@test IdealsQx[5]^2 == IdealsQx[5] * IdealsQx[5]
	@test ngens(IdealsQx[1]) == 3
	@test in(PolyQx[1], IdealsQx[2] : IdealsQx[3])
	@test dim(IdealsQx[1]) == 0
	@test dim(IdealsQx[1]) == 0

	@test issubset(intersect(IdealsNFx[3],IdealsNFx[4] + IdealsNFx[3]),(IdealsNFx[1]+IdealsNFx[2])*IdealsNFx[3])
	@test IdealsNFx[5]^2 == IdealsNFx[5] * IdealsNFx[5]
	@test ngens(IdealsNFx[1]) == 4
	@test in(PolyNFx[1], IdealsNFx[2] : IdealsNFx[3])
	@test dim(IdealsNFx[1]) == -1

	@test issubset(intersect(IdealsGFx[3],IdealsGFx[4] + IdealsGFx[3]),(IdealsGFx[1]+IdealsGFx[2])*IdealsGFx[3])
	@test IdealsGFx[4]^2 == IdealsGFx[4] * IdealsGFx[4]
	@test ngens(IdealsGFx[1]) == length(gens(IdealsGFx[1]))
	@test in(PolyGFx[1], IdealsGFx[2] : IdealsGFx[3])

	RFx,(x,y,z)= PolynomialRing(QQ,:x => 1:3)
	S, (a,b,c)=Singular.PolynomialRing(QQ, ["x$i" for i in 1:3])

	Oscar.getindex(IdealsQx[1].gens, Val(:S), 3)
	Oscar.iterate(IdealsQx[1].gens,8)
	Oscar.iterate(IdealsQx[1].gens,1)

	Singular.Fp(7)(Nemo.ResidueRing(ZZ,13)(3))
	Nemo.ResidueRing(ZZ,4)(Singular.Fp(7)(3))

	Oscar.singular_ring(RNmodx, :lex)
	Oscar.singular_ring(RNmodx, keep_ordering = false)
	Oscar.singular_ring(RNmodx, keep_ordering = true)
	Oscar.singular_ring(NFx[1], keep_ordering = true)
	Oscar.singular_ring(RFx)
	Oscar.singular_ring(RFx, keep_ordering = false)

	A=[]
	Oscar.push!(A, PolyQx[1])
	ideal(PolyNFx)
	ideal(A)
	ideal(parent(A[1]), A)
	@test gen(IdealsQx[1],1)==Poly[1]

	Oscar.oscar_assure(Oscar.MPolyIdeal(IdealsQx[1].gens.Ox, IdealsQx[1].gens.S))

	Oscar.groebner_basis(IdealsQx[1].gens)
	Oscar.groebner_basis(Oscar.BiPolyArray([x]))

	Oscar.groebner_basis(IdealsGFx[2].gens, ord = :lex, complete_reduction = false)
	####Oscar.syzygy_module(Poly)

	v=Singular.vector(S,convert(S,PolyQx[1]), convert(S,PolyQx[1]), convert(S, PolyQx[1]))
	Oscar.convert(Generic.FreeModule(S,3),v)
	Oscar.syzygy_generators([x  for x=PolyQx])

	Oscar.im_func(PolyQx[1], RFx, [2,3,1])
	T=Oscar.MPolyHom_vars{FmpqMPolyRing, FmpqMPolyRing}(RFx, RFx, [1,2,3])
	Oscar.MPolyHom_vars{FmpqMPolyRing, FmpqMPolyRing}(RFx, RFx, type = :names)
	T(PolyQx[1])
	Hecke.hom(RFx, RFx, [1,2,3])
	##Oscar._lift(IdealsNFx[1],IdealsNFx[1])
	Oscar.coordinates(PolyQx,PolyQx[1])

	##Oscar.eliminate(IdealsQx[1], [x,y])
	Oscar.eliminate(IdealsQx[1], [1,2,3])

	base_ring(IdealsQx[2])
	Oscar.groebner_basis(IdealsQx[2])
	Oscar.groebner_basis(IdealsQx[2], :lex)

	Q=Oscar.MPolyQuo(Qx,IdealsQx[1])
	gens(Q)
	ngens(Q)
	gen(Q,2)
	getindex(Q,2)
	f=Oscar.MPolyQuoElem{fmpq_mpoly}(PolyQx[1],Q)
	g=Oscar.MPolyQuoElem{fmpq_mpoly}(gens(IdealsQx[2])[1],Q)
	h=Oscar.MPolyQuoElem{fmpq_mpoly}(gens(IdealsQx[3])[1],Q)

	@test f*g-h*g == (f+ (-h))*g

	Oscar.mul!(g,h,f)
	Oscar.addeq!(f,g)
	Oscar.simplify!(f)
	@test !(f==g)
	Oscar.quo(Qx,IdealsQx[2])
	Oscar.lift(f)
	Q()
	Q(f)
	Q(Poly[1])
	zero(Q)
	Oscar._kbase(Q)
	Oscar.vector_space(QQ,Q)
end
