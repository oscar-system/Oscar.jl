using Test

@testset "ActionPolyRing - all tests" verbose = true begin
   
  __upr = Oscar.__upr
  __jtv = Oscar.__jtv
  __jtu_idx = Oscar.__jtu_idx
  __vtj = Oscar.__vtj
  __are_perms_up_to_date = Oscar.__are_perms_up_to_date
  __perm_for_sort = Oscar.__perm_for_sort
  __perm_for_sort_poly = Oscar.__perm_for_sort_poly
  __is_valid_jet = Oscar.__is_valid_jet
  __is_valid_partition = Oscar.__is_valid_partition
  
  @testset "Conformance tests - Ring interface" begin
    
    @testset "DifferencePolyRing - conformance tests" begin
      ConformanceTests.test_Ring_interface(difference_polynomial_ring(residue_ring(ZZ, ZZ(18))[1], 2, 3)[1])
      ConformanceTests.test_Ring_interface(difference_polynomial_ring(ZZ, 5, 8)[1])
    end
    
    @testset "DifferentialPolyRing - conformance tests" begin
      ConformanceTests.test_Ring_interface(differential_polynomial_ring(residue_ring(ZZ, ZZ(18))[1], 2, 3)[1])
      ConformanceTests.test_Ring_interface(differential_polynomial_ring(ZZ, 5, 8)[1])
    end

  end
  
  @testset "Check types" begin
    
    @test base_ring_type(DifferencePolyRing{QQFieldElem}) == QQField
    @test elem_type(DifferencePolyRing{QQFieldElem}) == DifferencePolyRingElem{QQFieldElem}
    @test parent_type(DifferencePolyRingElem{QQFieldElem}) == DifferencePolyRing{QQFieldElem}
    
    @test base_ring_type(DifferentialPolyRing{QQFieldElem}) == QQField
    @test elem_type(DifferentialPolyRing{QQFieldElem}) == DifferentialPolyRingElem{QQFieldElem}
    @test parent_type(DifferentialPolyRingElem{QQFieldElem}) == DifferentialPolyRing{QQFieldElem}

  end

  @testset "Check content" begin
    
    @testset "Failing constructions" begin
      for nelemsym in [-5, 0]
        @test_throws ArgumentError difference_polynomial_ring(ZZ, nelemsym, 1)
        @test_throws ArgumentError differential_polynomial_ring(ZZ, nelemsym, 1)
      end
      for ndiff in [-5, 0]
        @test_throws ArgumentError difference_polynomial_ring(ZZ, 1, ndiff)
        @test_throws ArgumentError differential_polynomial_ring(ZZ, 1, ndiff)
      end
      @test_throws ArgumentError difference_polynomial_ring(ZZ, Symbol[], 1)
      @test_throws ArgumentError differential_polynomial_ring(ZZ, Symbol[], 1)
    end

    @testset "Constructions with keywords" begin
      dpr0, vars0 = difference_polynomial_ring(ZZ, 3, 3; index_ordering_name = :invlex, partition_name = :pot)
      dpr1, vars1 = difference_polynomial_ring(ZZ, [:u1, :u2, :u3], 3; index_ordering_name = :invlex, partition_name = :pot)
      dpr2, vars2 = differential_polynomial_ring(ZZ, 3, 3; index_ordering_name = :invlex, partition_name = :pot)
      dpr3, vars3 = differential_polynomial_ring(ZZ, [:u1, :u2, :u3], 3; index_ordering_name = :invlex, partition_name = :pot)
      
      @testset "Check types" begin
        @test dpr0 isa DifferencePolyRing
        @test dpr1 isa DifferencePolyRing
        @test dpr2 isa DifferentialPolyRing
        @test dpr3 isa DifferentialPolyRing
        
        @test all(var -> var isa DifferencePolyRingElem, vars0)
        @test all(var -> var isa DifferencePolyRingElem, vars1)
        @test all(var -> var isa DifferentialPolyRingElem, vars2)
        @test all(var -> var isa DifferentialPolyRingElem, vars3)
      end

      for (dpr, vars) in [(dpr0, vars0), (dpr1, vars1), (dpr2, vars2), (dpr3, vars3)]
        u1, u2, u3 = vars[1], vars[2], vars[3]

        @testset "Internal getters" begin
          @test data(u1) === u1.p
          @test data(u2) === u2.p
          @test data(u3) === u3.p
          @test __upr(dpr) === dpr.upoly_ring
          @test __jtv(dpr) === dpr.jet_to_var
          @test __jtu_idx(dpr) === dpr.jet_to_upoly_idx
          @test __are_perms_up_to_date(dpr) === dpr.are_perms_up_to_date
          @test __perm_for_sort(dpr) === dpr.permutation
          @test __perm_for_sort_poly(u1) === u1.permutation
          @test __perm_for_sort_poly(u2) === u2.permutation
        end

        @testset "Check internals at construction" begin
          upr = __upr(dpr)
          u1p, u2p, u3p = data(u1), data(u2), data(u3)
          @test parent(u1p) === upr
          @test parent(u2p) === upr
          @test parent(u3p) === upr
          @test gens(upr) == [u1p, u2p, u3p]
          
          #Dictionaries
          jtu = __jtu_idx(dpr)
          jtv = __jtv(dpr)
          vtj = __vtj(dpr)
          for key in keys(jtv)
            @test vtj[jtv[key]] === key
          end
          for key in keys(vtj)
            @test jtv[vtj[key]] === key
          end
          @test Set(keys(jtv)) == Set([(1, [0,0,0]), (2, [0,0,0]), (3, [0,0,0])])
          @test Set(keys(vtj)) == Set([u1, u2, u3])
          @test keys(jtu) == keys(jtv)
          @test jtu[(1, [0,0,0])] == 1 && jtu[(2, [0,0,0])] == 2 && jtu[(3, [0,0,0])] == 3
          
          @test u1 == jtv[(1, [0,0,0])]
          @test u2 == jtv[(2, [0,0,0])]
          @test u3 == jtv[(3, [0,0,0])]

          @test u1 !== jtv[(1, [0,0,0])]
          @test u2 !== jtv[(2, [0,0,0])]
          @test u3 !== jtv[(3, [0,0,0])]

          #Permutations
          @test __perm_for_sort_poly(u1) == [1]
          @test __perm_for_sort_poly(u2) == [1]
          @test __perm_for_sort_poly(u3) == [1]

          @test u1 > u2 > u3
          @test __perm_for_sort(dpr) == [1,2,3]
        end
        
        @testset "Check other" begin
          @test !is_univariate(dpr)
          @test all(var -> is_univariate(var), vars)
          @test_throws ErrorException to_univariate(u1*u2)
          @test typeof(to_univariate(u1)) == ZZPolyRingElem
          R, univ_var = polynomial_ring(ZZ; cached = false)
          @test to_univariate(R, u1^2 - u1) == univ_var^2 - univ_var
          @test parent(to_univariate(R, u1^2 - u1)) === R
          
          @test all(var -> is_irreducible(var), vars)
          @test !is_irreducible(u1^2)
          @test !is_irreducible(u1*u2)
          @test !is_irreducible(zero(dpr))
          @test !is_irreducible(one(dpr))
          
          @test all(var -> is_monomial(var), vars)
          @test all(var -> is_term(var), vars)
          @test all(var -> !is_monomial(3*var), vars)
          @test all(var -> is_term(3*var), vars)
          @test all(var -> !is_monomial(3*var+1), vars)
          @test all(var -> !is_term(3*var+1), vars)
          
          f = 2 * (u1 - 1) * (u1 + 1)
          g = -5*u1^3 * u2^2
          factf, factg = factor(f), factor(g)
          factsf, factsg = factor_squarefree(f), factor_squarefree(g)
          @test unit(factf) == ZZ(1)
          @test unit(factsf) == ZZ(1)
          @test Set(factf) == Set([Pair(dpr(2), 1), Pair(u1 - 1, 1), Pair(u1 + 1, 1)])
          @test Set(factsf) == Set([Pair(dpr(2), 1), Pair(u1^2 - 1, 1)])
          @test unit(factg) == ZZ(-1)
          @test unit(factsg) == ZZ(-1)
          @test Set(factg) == Set([Pair(dpr(5), 1), Pair(u1, 3), Pair(u2, 2)])
          @test Set(factsg) == Set([Pair(dpr(5), 1), Pair(u1, 3), Pair(u2, 2)])
        end
        
        @testset "Check public fields at construction" begin
          @test base_ring(dpr) == ZZ
          @test all(var -> base_ring(var) == ZZ, vars)
          @test elementary_symbols(dpr) == [:u1, :u2, :u3]
          @test ndiffs(dpr) == 3
          @test n_elementary_symbols(dpr) == 3
          @test all(var -> parent(var) === dpr, vars)

          ran = ranking(dpr)
          if dpr isa DifferencePolyRing
            @test ran isa ActionPolyRingRanking{DifferencePolyRing{ZZRingElem}}
          end
          if dpr isa DifferentialPolyRing
            @test ran isa ActionPolyRingRanking{DifferentialPolyRing{ZZRingElem}}
          end
          @test partition(ran) == [[1,0,0], [0,1,0], [0,0,1]]
          @test index_ordering_matrix(ran) == ZZ[0 0 1; 0 1 0; 1 0 0]
          @test riquier_matrix(ran) == ZZ[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0; 0 0 0 1 0 0]
        end
        
        # Let's add some variables
        @testset "Check adding variables" begin
          for i in [1,2,3]
            @test __is_valid_jet(dpr, i, [17,0,3])
          end
          @test !__is_valid_jet(dpr, 0, [17,0,3])
          @test !__is_valid_jet(dpr, 4, [17,0,3])

          for v in [[11,0,5], [2,7,1], [4,1,4]]
            @test __is_valid_jet(dpr, 1, v)
          end
          @test !__is_valid_jet(dpr, 1, [1,1])
          @test !__is_valid_jet(dpr, 1, [1,1,1,1])
          @test !__is_valid_jet(dpr, 1, [0,0,-5])

          @test number_of_generators(dpr) == 3
          @test ngens(dpr) == 3
          @test number_of_variables(dpr) == 3
          @test nvars(dpr) == 3
          @test_throws ArgumentError dpr[0,[0,0,0]]
          @test_throws ArgumentError dpr[4,[0,0,0]]
          @test_throws ArgumentError dpr[1,[10,4,-1]]
          u1_100 = dpr[1, [1,0,0]]
          @test ngens(dpr) == 4
          @test nvars(dpr) == 4
          u1_010 = dpr[(1, [0,1,0])]
          @test ngens(dpr) == 5
          @test nvars(dpr) == 5
          u2_100 = gen(dpr, 2, [1,0,0])
          @test ngens(dpr) == 6
          @test nvars(dpr) == 6
          gen(dpr, (2, [1,0,0])) #Nothing new gets added here
          @test ngens(dpr) == 6
          @test nvars(dpr) == 6
          new_vars_idx = [(1, [1,0,0]), (1, [0,1,0]), (3, [1,1,1]), (2, [1,0,0])]
          gens(dpr, new_vars_idx)
          u3_111 = dpr[3, [1,1,1]]
          @test ngens(dpr) == 7
          @test nvars(dpr) == 7

          @test u1_010 > u1_100 > u1 > u2_100 > u2 > u3_111 > u3 #position over term and invlex
          @test gens(dpr) == [u1_010, u1_100, u1, u2_100, u2, u3_111, u3]
          @test all(var -> is_gen(var), gens(dpr))
            
          for i in 1:ngens(dpr)
            @test gen(dpr, i) == gens(dpr)[i]
            @test __jtv(dpr)[__vtj(dpr)[gen(dpr, i)]] !== gen(dpr, i)
            Rtmp = polynomial_ring(ZZ; cached = false)[1]
            @test to_univariate(Rtmp, gen(dpr, i)) == gen(Rtmp)
            @test to_univariate(gen(dpr, i)) isa ZZPolyRingElem
          end
        
          @testset "Check internals after adding variables" begin
            upr = __upr(dpr)
            @test all(var -> parent(data(var)) === upr, gens(dpr))
          
            #Dictionaries
            jtu = __jtu_idx(dpr)
            jtv = __jtv(dpr) 
            vtj = __vtj(dpr)

            for key in keys(jtv)
              @test vtj[jtv[key]] === key
            end
            for key in keys(vtj)
              @test jtv[vtj[key]] === key
            end
            @test Set(keys(jtv)) == Set([(1, [0,0,0]), (2, [0,0,0]), (3, [0,0,0]), (1, [0,1,0]), (1, [1,0,0]), (2, [1,0,0]), (3,[1,1,1])])
            @test Set(keys(vtj)) == Set([u1, u2, u3, u1_010, u1_100, u2_100, u3_111])
            @test keys(jtu) == keys(jtv)
            @test jtu[(1, [0,0,0])] == 1 && jtu[(2, [0,0,0])] == 2 && jtu[(3, [0,0,0])] == 3 && jtu[(1, [0,1,0])] == 5 && jtu[(1, [1,0,0])] == 4 && jtu[(2, [1,0,0])] == 6 && jtu[(3, [1,1,1])] == 7
          
            #Permutations
            @test !__are_perms_up_to_date(dpr)
            @test __perm_for_sort(dpr) == [5,4,1,6,2,7,3]
            @test __are_perms_up_to_date(dpr)
          end
        
          @testset "Check public fields after adding variables" begin
            @test base_ring(dpr) == ZZ
            @test all(var -> base_ring(var) == ZZ, vars)
            @test elementary_symbols(dpr) == [:u1, :u2, :u3]
            @test ndiffs(dpr) == 3
            @test n_elementary_symbols(dpr) == 3

            ran = ranking(dpr)
            if dpr isa DifferencePolyRing
              @test ran isa ActionPolyRingRanking{DifferencePolyRing{ZZRingElem}}
            end
            if dpr isa DifferentialPolyRing
              @test ran isa ActionPolyRingRanking{DifferentialPolyRing{ZZRingElem}}
            end
            @test partition(ran) == [[1,0,0], [0,1,0], [0,0,1]]
            @test index_ordering_matrix(ran) == ZZ[0 0 1; 0 1 0; 1 0 0]
            @test riquier_matrix(ran) == ZZ[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0; 0 0 0 1 0 0]
          end

          @testset "Constant polynomials" begin
            @test dpr() == zero(dpr)
            @test dpr(0) == zero(dpr)
            @test is_zero(dpr())
            @test !is_unit(dpr())
            @test is_constant(dpr())
            @test is_univariate(dpr())
            @test_throws ErrorException inv(dpr())
            @test_throws ArgumentError trailing_coefficient(dpr())
            @test_throws ArgumentError trailing_term(dpr())
            @test_throws ArgumentError trailing_monomial(dpr())
            @test_throws ArgumentError leading_coefficient(dpr())
            @test_throws ArgumentError leading_term(dpr())
            @test_throws ArgumentError leading_monomial(dpr())
            @test length(dpr()) == 0
            @test collect(coefficients(dpr())) == ZZRingElem[]
            @test collect(exponents(dpr())) == Vector{Int}[]
            @test collect(monomials(dpr())) == elem_type(dpr)[]
            @test collect(terms(dpr())) == elem_type(dpr)[]
            @test degrees(dpr()) == fill(-1, nvars(dpr))
            for i in 1:length(nvars(dpr))
              @test degree(dpr(), i) == degrees(dpr())[i]
            end
            @test_throws BoundsError degree(dpr(), 0)
            @test_throws BoundsError degree(dpr(), nvars(dpr) + 1)
            @test total_degree(dpr()) == -1

            @test dpr(1) == one(dpr)
            @test is_one(dpr(1))
            @test is_unit(dpr(1))
            @test is_constant(dpr(1))
            @test is_univariate(dpr(1))
            @test inv(dpr(1)) == 1
            @test trailing_coefficient(dpr(1)) == one(ZZ)
            @test trailing_term(dpr(1)) == dpr(1)
            @test trailing_monomial(dpr(1)) == dpr(1)
            @test leading_coefficient(dpr(1)) == one(ZZ)
            @test leading_term(dpr(1)) == dpr(1)
            @test leading_monomial(dpr(1)) == dpr(1)
            @test length(dpr(1)) == 1
            @test collect(coefficients(dpr(1))) == ZZRingElem[1]
            @test collect(exponents(dpr(1))) == Vector{Int}[fill(0, nvars(dpr))]
            @test collect(monomials(dpr(1))) == elem_type(dpr)[dpr(1)]
            @test collect(terms(dpr(1))) == elem_type(dpr)[dpr(1)]
            @test degrees(dpr(1)) == fill(0, nvars(dpr))
            for i in 1:nvars(dpr)
              @test degree(dpr(1), i) == degrees(dpr(1))[i]
            end
            @test_throws BoundsError degree(dpr(1), 0)
            @test_throws BoundsError degree(dpr(1), nvars(dpr) + 1)
            @test total_degree(dpr(1)) == 0 

            @test !is_unit(dpr(-2))
            @test is_constant(dpr(-2))
            @test is_univariate(dpr(-2))
            @test_throws ErrorException inv(dpr(-2))
            @test trailing_coefficient(dpr(-2)) == ZZ(-2)
            @test trailing_term(dpr(-2)) == dpr(-2)
            @test trailing_monomial(dpr(-2)) == dpr(1)
            @test leading_coefficient(dpr(-2)) == ZZ(-2)
            @test leading_term(dpr(-2)) == dpr(-2)
            @test leading_monomial(dpr(-2)) == dpr(1)
            @test length(dpr(-2)) == 1
            @test collect(coefficients(dpr(-2))) == ZZRingElem[-2]
            @test collect(exponents(dpr(-2))) == Vector{Int}[fill(0, nvars(dpr))]
            @test collect(monomials(dpr(-2))) == elem_type(dpr)[dpr(1)]
            @test collect(terms(dpr(-2))) == elem_type(dpr)[dpr(-2)]
            @test degrees(dpr(-2)) == fill(0, nvars(dpr))
            for i in 1:nvars(dpr)
              @test degree(dpr(-2), i) == degrees(dpr(-2))[i]
            end
            @test_throws BoundsError degree(dpr(-2), 0)
            @test_throws BoundsError degree(dpr(-2), nvars(dpr) + 1)
            @test total_degree(dpr(-2)) == 0 

            @test dpr(-5) == -dpr(5)
            @test -1 * dpr(5) == dpr(-5)
            @test ZZ(-1) * dpr(5) == dpr(-5)
            @test dpr(5) + dpr(-5) == dpr(5 + -5)
            @test dpr(5)^2 == dpr(5) * dpr(5)
            @test dpr(25) == dpr(5) * dpr(5)
            @test is_square_with_sqrt(dpr(25)) == (true, 5)
            @test !is_square(dpr(5))
            @test divexact(dpr(5), -5) == -1
            @test divexact(dpr(5), ZZ(-5)) == -1
            @test divexact(dpr(5), dpr(-5)) == -1
            @test_throws ErrorException divexact(dpr(5), 2)
            @test_throws ErrorException divexact(dpr(5), ZZ(2))
            @test_throws ErrorException divexact(dpr(5), dpr(2))
          end

          @testset "Non-constant polynomials" begin
            #Recall that u1_010 > u1_100 > u1 > u2_100 > u2 > u3_111 > u3:
            @test u1_010 > u1_100 > u1 > u2_100 > u2 > u3_111 > u3 #position over term and invlex
            @test [var_index(u1_010), var_index(u1), var_index(u3_111)] == [1,3,6]
            f = (u3_111 - 2 * u2_100) * (u1 - u1_100 + 3)
            g = (u1_010 - 2) * (u1 - u1_100 + 3)
            
            @test f == 2*u2_100*u1_100 - u3_111*u1_100 - 2*u2_100*u1 + u3_111*u1 - 6*u2_100 + 3*u3_111
            @test length(f) == 6
            @test __perm_for_sort_poly(f) == [3,4,1,2,5,6]

            cf = collect(coefficients(f))
            ef = collect(exponents(f))
            mf = collect(monomials(f))
            tf = collect(terms(f))
            @test typeof(cf) == Vector{ZZRingElem}
            @test typeof(ef) == Vector{Vector{Int}}
            @test typeof(mf) == Vector{typeof(f)}
            @test typeof(tf) == Vector{typeof(f)}
            
            @test cf == ZZRingElem[2,-1,-2,1,-6,3]
            @test ef == [[0,1,0,1,0,0,0],[0,1,0,0,0,1,0],[0,0,1,1,0,0,0],[0,0,1,0,0,1,0],[0,0,0,1,0,0,0],[0,0,0,0,0,1,0]]
            @test mf == [u2_100*u1_100, u3_111*u1_100, u2_100*u1, u3_111*u1, u2_100, u3_111]
            @test tf == [2*u2_100*u1_100, -u3_111*u1_100, -2*u2_100*u1, u3_111*u1, -6*u2_100, 3*u3_111]
          
            @test !is_zero(f)
            @test !is_unit(f)
            @test !is_constant(dpr(f))
            @test !is_univariate(f)
            @test_throws ErrorException to_univariate(f)
            @test_throws ErrorException to_univariate(polynomial_ring(ZZ; cached = false)[1], f)
            @test_throws ErrorException inv(f)
            @test trailing_coefficient(f) == ZZ(3)
            @test trailing_monomial(f) == u3_111
            @test trailing_term(f) == 3*u3_111 
            @test leading_coefficient(f) == ZZ(2)
            @test leading_monomial(f) == u2_100*u1_100
            @test leading_term(f) == 2*u2_100*u1_100
            @test tail(f) == f - leading_term(f)
            @test degrees(f) == [0,1,1,1,0,1,0]
            for i in 1:length(nvars(dpr))
              @test degree(f, i) == degrees(f)[i]
              @test degree(f, i) == degree(f, gen(parent(f), i))
            end
            @test degree(f, 3, [5,5,5]) == 0 && ngens(dpr) == 7
            @test_throws ArgumentError degree(f, 0, [5,5,5])
            @test_throws ArgumentError degree(f, 1, [5,5,5,5])
            @test_throws BoundsError degree(f, 0)
            @test_throws BoundsError degree(f, nvars(dpr) + 1)
            @test total_degree(f) == 2
            
            @test leader(f) == u1_100
            @test degree(f, leader(f)) == 1
            @test initial(f) == leading_term(f)/leader(f)

            @test dpr(data(f) + data(g)) == f + g
            @test dpr(-data(f)) == -f
            @test dpr(-1*data(f)) == -1 * f
            @test f^2 == f * f
            @test dpr(data(f)^2) == f^2
            @test dpr(data(f) * data(g)) == f * g
            @test_throws ErrorException f/g 
            @test_throws ErrorException divexact(f, g)
            @test is_one(f/f)
            @test is_one(divexact(f,f))

            #g = (u1_010 - 2) * (u1 - u1_100 + 3)
            @test g == -u1_010*u1_100 + u1_010*u1 + 3*u1_010 + 2*u1_100 - 2*u1 - 6
            @test length(g) == 6
            @test __perm_for_sort_poly(g) == [3,1,5,4,2,6]
            cg = collect(coefficients(g))
            eg = collect(exponents(g))
            mg = collect(monomials(g))
            tg = collect(terms(g))
           
            @test cg == ZZRingElem[-1,1,3,2,-2,-6]
            @test eg == [[1,1,0,0,0,0,0],[1,0,1,0,0,0,0],[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,0,0,0,0]]
            @test mg == [u1_010*u1_100, u1_010*u1, u1_010, u1_100, u1, dpr(1)]
            @test tg == [-1*u1_010*u1_100, u1_010*u1, 3*u1_010, 2*u1_100, -2*u1, -6*dpr(1)]
            
            @test leader(g) == u1_010
            @test degree(g, leader(g)) == 1
            @test initial(g) == leading_term(g)/leader(g)

          end
        
        end #End check adding variables

        @testset "Check change ranking" begin
          set_ranking!(dpr; partition = [[0,1,1], [1,0,0]], index_ordering_name = :degrevlex)
         
          @test ngens(dpr) == 7
          @test nvars(dpr) == 7
          u1_100 = dpr[1, [1,0,0]]
          u1_010 = dpr[(1, [0,1,0])]
          u2_100 = gen(dpr, 2, [1,0,0])
          u3_111 = dpr[3, [1,1,1]]
          @test ngens(dpr) == 7
          @test nvars(dpr) == 7
          
          # u1_010 > u1_100 > u1 > u2_100 > u2 > u3_111 > u3 #previously
          @test u3_111 > u2_100 > u2 > u3 > u1_100 > u1_010 > u1
          @test [var_index(u1_010), var_index(u1), var_index(u3_111)] == [6,7,1]       
          @test gens(dpr) == [u3_111, u2_100, u2, u3, u1_100, u1_010, u1]
          @test all(var -> is_gen(var), gens(dpr))
            
          for i in 1:ngens(dpr)
            @test gen(dpr, i) == gens(dpr)[i]
            @test __jtv(dpr)[__vtj(dpr)[gen(dpr, i)]] !== gen(dpr, i)
            Rtmp = polynomial_ring(ZZ; cached = false)[1]
            @test to_univariate(Rtmp, gen(dpr, i)) == gen(Rtmp)
            @test to_univariate(gen(dpr, i)) isa ZZPolyRingElem
          end
          
          @testset "Check internals after changing ranking" begin
            upr = __upr(dpr)
            @test all(var -> parent(data(var)) === upr, gens(dpr))
          
            #Dictionaries
            jtu = __jtu_idx(dpr)
            jtv = __jtv(dpr) 
            vtj = __vtj(dpr)

            for key in keys(jtv)
              @test vtj[jtv[key]] === key
            end
            for key in keys(vtj)
              @test jtv[vtj[key]] === key
            end
            @test Set(keys(jtv)) == Set([(1, [0,0,0]), (2, [0,0,0]), (3, [0,0,0]), (1, [0,1,0]), (1, [1,0,0]), (2, [1,0,0]), (3,[1,1,1])])
            @test Set(keys(vtj)) == Set([u1, u2, u3, u1_010, u1_100, u2_100, u3_111])
            @test keys(jtu) == keys(jtv)
            @test jtu[(1, [0,0,0])] == 1 && jtu[(2, [0,0,0])] == 2 && jtu[(3, [0,0,0])] == 3 && jtu[(1, [0,1,0])] == 5 && jtu[(1, [1,0,0])] == 4 && jtu[(2, [1,0,0])] == 6 && jtu[(3, [1,1,1])] == 7
          
            #Permutations
            @test !__are_perms_up_to_date(dpr)
            @test __perm_for_sort(dpr) == [7,6,2,3,4,5,1]
            @test __are_perms_up_to_date(dpr)
          end
        
          @testset "Check public fields after changing ranking" begin
            @test base_ring(dpr) == ZZ
            @test all(var -> base_ring(var) == ZZ, vars)
            @test elementary_symbols(dpr) == [:u1, :u2, :u3]
            @test ndiffs(dpr) == 3
            @test n_elementary_symbols(dpr) == 3
          
            ran = ranking(dpr)
            if dpr isa DifferencePolyRing
              @test ran isa ActionPolyRingRanking{DifferencePolyRing{ZZRingElem}}
            end
            if dpr isa DifferentialPolyRing
              @test ran isa ActionPolyRingRanking{DifferentialPolyRing{ZZRingElem}}
            end
            @test partition(ran) == [[0,1,1], [1,0,0]]
            @test index_ordering_matrix(ran) == ZZ[1 1 1; 0 0 -1; 0 -1 0]
            @test riquier_matrix(ran) == ZZ[0 1 1 0 0 0; 0 0 0 1 1 1; 0 0 0 0 0 -1; 0 0 0 0 -1 0; 0 1 0 0 0 0]
          end

          @testset "Constant polynomials" begin
            @test dpr() == zero(dpr)
            @test dpr(0) == zero(dpr)
            @test dpr() == ZZ()
            @test is_zero(dpr())
            @test !is_unit(dpr())
            @test is_constant(dpr())
            @test is_univariate(dpr())
            @test_throws ErrorException inv(dpr())
            @test_throws ArgumentError trailing_coefficient(dpr())
            @test_throws ArgumentError trailing_term(dpr())
            @test_throws ArgumentError trailing_monomial(dpr())
            @test_throws ArgumentError leading_coefficient(dpr())
            @test_throws ArgumentError leading_term(dpr())
            @test_throws ArgumentError leading_monomial(dpr())
            @test length(dpr()) == 0
            @test collect(coefficients(dpr())) == ZZRingElem[]
            @test collect(exponents(dpr())) == Vector{Int}[]
            @test collect(monomials(dpr())) == elem_type(dpr)[]
            @test collect(terms(dpr())) == elem_type(dpr)[]
            @test degrees(dpr()) == fill(-1, nvars(dpr))
            for i in 1:length(nvars(dpr))
              @test degree(dpr(), i) == degrees(dpr())[i]
            end
            @test_throws BoundsError degree(dpr(), 0)
            @test_throws BoundsError degree(dpr(), nvars(dpr) + 1)
            @test total_degree(dpr()) == -1
            @test_throws ArgumentError leader(dpr())
            @test initial(dpr()) == ZZ()

            @test dpr(1) == one(dpr)
            @test dpr(1) == ZZ(1)
            @test is_one(dpr(1))
            @test is_unit(dpr(1))
            @test is_constant(dpr(1))
            @test is_univariate(dpr(1))
            @test inv(dpr(1)) == 1
            @test trailing_coefficient(dpr(1)) == one(ZZ)
            @test leading_coefficient(dpr(1)) == one(ZZ)
            @test leading_term(dpr(1)) == dpr(1)
            @test leading_monomial(dpr(1)) == dpr(1)
            @test length(dpr(1)) == 1
            @test collect(coefficients(dpr(1))) == ZZRingElem[1]
            @test collect(exponents(dpr(1))) == Vector{Int}[fill(0, nvars(dpr))]
            @test collect(monomials(dpr(1))) == elem_type(dpr)[dpr(1)]
            @test collect(terms(dpr(1))) == elem_type(dpr)[dpr(1)]
            @test degrees(dpr(1)) == fill(0, nvars(dpr))
            for i in 1:nvars(dpr)
              @test degree(dpr(1), i) == degrees(dpr(1))[i]
            end
            @test_throws BoundsError degree(dpr(1), 0)
            @test_throws BoundsError degree(dpr(1), nvars(dpr) + 1)
            @test total_degree(dpr(1)) == 0 
            @test_throws ArgumentError leader(dpr(1))
            @test initial(dpr(1)) == ZZ(1)

            @test dpr(-2) == ZZ(-2)
            @test !is_unit(dpr(-2))
            @test is_constant(dpr(-2))
            @test is_univariate(dpr(-2))
            @test_throws ErrorException inv(dpr(-2))
            @test trailing_coefficient(dpr(-2)) == ZZ(-2)
            @test leading_coefficient(dpr(-2)) == ZZ(-2)
            @test leading_term(dpr(-2)) == dpr(-2)
            @test leading_monomial(dpr(-2)) == dpr(1)
            @test length(dpr(-2)) == 1
            @test collect(coefficients(dpr(-2))) == ZZRingElem[-2]
            @test collect(exponents(dpr(-2))) == Vector{Int}[fill(0, nvars(dpr))]
            @test collect(monomials(dpr(-2))) == elem_type(dpr)[dpr(1)]
            @test collect(terms(dpr(-2))) == elem_type(dpr)[dpr(-2)]
            @test degrees(dpr(-2)) == fill(0, nvars(dpr))
            for i in 1:nvars(dpr)
              @test degree(dpr(-2), i) == degrees(dpr(-2))[i]
            end
            @test_throws BoundsError degree(dpr(-2), 0)
            @test_throws BoundsError degree(dpr(-2), nvars(dpr) + 1)
            @test total_degree(dpr(-2)) == 0 
            @test_throws ArgumentError leader(dpr(-2))
            @test initial(dpr(-2)) == ZZ(-2)

            @test dpr(-5) == -dpr(5)
            @test -1 * dpr(5) == dpr(-5)
            @test ZZ(-1) * dpr(5) == dpr(-5)
            @test dpr(5) + dpr(-5) == dpr(5 + -5)
            @test dpr(5)^2 == dpr(5) * dpr(5)
            @test dpr(25) == dpr(5) * dpr(5)
            @test is_square_with_sqrt(dpr(25)) == (true, 5)
            @test !is_square(dpr(5))
            @test divexact(dpr(5), -5) == -1
            @test divexact(dpr(5), ZZ(-5)) == -1
            @test divexact(dpr(5), dpr(-5)) == -1
            @test_throws ErrorException divexact(dpr(5), 2)
            @test_throws ErrorException divexact(dpr(5), ZZ(2))
            @test_throws ErrorException divexact(dpr(5), dpr(2))
          end
          
          @testset "Non-constant polynomials" begin
            #Recall that u3_111 > u2_100 > u2 > u3 > u1_100 > u1_010 > u1
            f = (u3_111 - 2 * u2_100) * (u1 - u1_100 + 3)
            g = (u1_010 - 2) * (u1 - u1_100 + 3)
            @test dpr(data(f)) == f
            @test dpr(data(f)) !== f
            @test dpr(data(data(f))) == f
            @test dpr(data(data(f))) !== f
            @test dpr(f) === f

            @test dpr(data(g)) == g
            @test dpr(data(g)) !== g
            @test dpr(data(data(g))) == g
            @test dpr(data(data(g))) !== g
            @test dpr(g) === g

            @test f == 2*u2_100*u1_100 - u3_111*u1_100 - 2*u2_100*u1 + u3_111*u1 - 6*u2_100 + 3*u3_111
            @test length(f) == 6
            @test __perm_for_sort_poly(f) == [4,2,6,3,1,5]
            cf = collect(coefficients(f))
            ef = collect(exponents(f))
            mf = collect(monomials(f))
            tf = collect(terms(f))
            @test typeof(cf) == Vector{ZZRingElem}
            @test typeof(ef) == Vector{Vector{Int}}
            @test typeof(mf) == Vector{typeof(f)}
            @test typeof(tf) == Vector{typeof(f)}
            
            @test cf == ZZRingElem[-1,1,3,2,-2,-6]
            @test ef == [[1,0,0,0,1,0,0],[1,0,0,0,0,0,1],[1,0,0,0,0,0,0],[0,1,0,0,1,0,0],[0,1,0,0,0,0,1],[0,1,0,0,0,0,0]]
            @test mf == [u3_111*u1_100, u3_111*u1, u3_111, u2_100*u1_100, u2_100*u1, u2_100]
            @test tf == [-u3_111*u1_100, u3_111*u1, 3*u3_111, 2*u2_100*u1_100, -2*u2_100*u1, -6*u2_100]
            
            @test !is_zero(f)
            @test !is_unit(f)
            @test !is_constant(dpr(f))
            @test !is_univariate(f)
            @test_throws ErrorException to_univariate(f)
            @test_throws ErrorException to_univariate(polynomial_ring(ZZ; cached = false)[1], f)
            @test_throws ErrorException inv(f)
            @test trailing_coefficient(f) == ZZ(-6)
            @test trailing_monomial(f) == u2_100
            @test trailing_term(f) == -6*u2_100
            @test leading_coefficient(f) == ZZ(-1)
            @test leading_monomial(f) == u3_111*u1_100
            @test leading_term(f) == -1*u3_111*u1_100
            @test tail(f) == f - leading_term(f)
            @test degrees(f) == [1,1,0,0,1,0,1]
            for i in 1:length(nvars(dpr))
              @test degree(f, i) == degrees(f)[i]
              @test degree(f, i) == degree(f, gen(parent(f), i))
            end
            @test_throws BoundsError degree(f, 0)
            @test_throws BoundsError degree(f, nvars(dpr) + 1)
            @test total_degree(f) == 2
            
            @test leader(f) == u3_111
            @test degree(f, leader(f)) == 1
            @test initial(f) == leading_term(f)/leader(f)

            @test dpr(data(f) + data(g)) == f + g
            @test dpr(-data(f)) == -f
            @test dpr(-1*data(f)) == -1 * f
            @test f^2 == f * f
            @test dpr(data(f)^2) == f^2
            @test dpr(data(f) * data(g)) == f * g
            @test_throws ErrorException f/g 
            @test_throws ErrorException divexact(f, g)
            @test is_one(f/f)
            @test is_one(divexact(f,f))

            #g = (u1_010 - 2) * (u1 - u1_100 + 3)
            @test g == -u1_010*u1_100 + u1_010*u1 + 3*u1_010 + 2*u1_100 - 2*u1 - 6
            @test length(g) == 6
            @test __perm_for_sort_poly(g) == [3,4,1,5,2,6]
            cg = collect(coefficients(g))
            eg = collect(exponents(g))
            mg = collect(monomials(g))
            tg = collect(terms(g))
            
            @test leader(g) == u1_100
            @test degree(g, leader(g)) == 1
            @test initial(g) == leading_term(g)/leader(g)
            
            @test cg == ZZRingElem[-1,2,1,3,-2,-6]
            @test eg == [[0,0,0,0,1,1,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,1],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1],[0,0,0,0,0,0,0]]
            @test mg == [u1_010*u1_100, u1_100, u1_010*u1, u1_010, u1, dpr(1)]
            @test tg == [-u1_010*u1_100, 2*u1_100, u1_010*u1, 3*u1_010, -2*u1, -6*dpr(1)]

          end

        end #End check change ranking

      end #End for loop
      
    end #Construction with keywords

    @testset "Further construction" begin
      dpr0, vars0 = difference_polynomial_ring(QQ, 3, 3)
      dpr1, vars1 = difference_polynomial_ring(QQ, [:u1, :u2, :u3], 3)
      dpr2, vars2 = differential_polynomial_ring(QQ, 3, 3)
      dpr3, vars3 = differential_polynomial_ring(QQ, [:u1, :u2, :u3], 3)
      
      for (dpr, vars) in [(dpr0, vars0), (dpr1, vars1), (dpr2, vars2), (dpr3, vars3)]
        u1, u2, u3 = vars[1], vars[2], vars[3]

        f = u1*u2
        g = -3*u1^2*u3 + 4*u2
        @test gens(dpr) == [u1, u2, u3] #order of variables

        @testset "derivative" begin
          for i in 1:ngens(dpr)
            @test is_zero(derivative(dpr(), i))
            @test is_zero(derivative(dpr(1), i))
            @test is_zero(derivative(dpr(-2), i))
            @test derivative(g, i) == derivative(g, gen(dpr, i))
            @test derivative(g, i) == derivative(g, (i, [0,0,0]))
          end
          @test derivative(g, u1) == -6*u1*u3
          @test derivative(g, u2) == dpr(4)
          @test derivative(g, u3) == -3*u1^2
          @test_throws BoundsError derivative(g, 0)
          @test_throws BoundsError derivative(g, ngens(dpr) + 1)
          @test_throws ArgumentError derivative(g, (0, [1,1,1]))
          @test_throws ArgumentError derivative(g, (1, [1,1,1,1]))
          @test_throws ArgumentError derivative(g, (1, [1,1]))
          @test_throws ArgumentError derivative(g, (1, [1,-1,1]))
          @test is_zero(derivative(g, (2, [10,4,2]))) && ngens(dpr) == 3
        end

        @testset "diff action" begin
          if dpr isa DifferencePolyRing
            @test is_zero(diff_action(dpr(), 1))
            @test is_zero(diff_action(dpr(), ndiffs(dpr)))
            @test_throws ArgumentError diff_action(dpr(), 0)
            @test_throws ArgumentError diff_action(dpr(), ndiffs(dpr) + 1)
            @test diff_action(dpr(-2), 1) == dpr(-2)
            @test diff_action(dpr(-2), [0,0,0]) == dpr(-2)
            @test_throws ArgumentError diff_action(dpr(-2), [1,1,1,1]) 
            @test_throws ArgumentError diff_action(dpr(-2), [1,1]) 
            @test_throws ArgumentError diff_action(dpr(-2), [1,-1,1]) 
            
            @test ngens(dpr) == 3
            @test diff_action(f, 1) == dpr[1, [1,0,0]] * dpr[2, [1,0,0]]
            @test ngens(dpr) == 5
            @test diff_action(f, 2) == dpr[1, [0,1,0]] * dpr[2, [0,1,0]]
            @test ngens(dpr) == 7
            @test diff_action(f, 3) == dpr[1, [0,0,1]] * dpr[2, [0,0,1]]
            @test ngens(dpr) == 9

            @test diff_action(g, 1) == -3*dpr[1, [1,0,0]]^2 * dpr[3, [1,0,0]] + 4*dpr[2, [1,0,0]]
            @test ngens(dpr) == 10
            @test diff_action(g, 2) == -3*dpr[1, [0,1,0]]^2 * dpr[3, [0,1,0]] + 4*dpr[2, [0,1,0]]
            @test ngens(dpr) == 11
            @test diff_action(g, 3) == -3*dpr[1, [0,0,1]]^2 * dpr[3, [0,0,1]] + 4*dpr[2, [0,0,1]]
            @test ngens(dpr) == 12

            @test diff_action(f, [4,5,6]) == dpr[1, [4,5,6]] * dpr[2, [4,5,6]]
            @test ngens(dpr) == 14
            @test diff_action(g, [4,5,6]) == -3*dpr[1, [4,5,6]]^2 * dpr[3, [4,5,6]] + 4*dpr[2, [4,5,6]]
            @test ngens(dpr) == 15
          end
          if dpr isa DifferentialPolyRing
            @test is_zero(diff_action(dpr(), 1))
            @test is_zero(diff_action(dpr(), ndiffs(dpr)))
            @test_throws ArgumentError diff_action(dpr(), 0)
            @test_throws ArgumentError diff_action(dpr(), ndiffs(dpr) + 1)
            @test is_zero(diff_action(dpr(-2), 1))
            @test diff_action(dpr(-2), [0,0,0]) == -2
            @test_throws ArgumentError diff_action(dpr(-2), [1,1,1,1]) 
            @test_throws ArgumentError diff_action(dpr(-2), [1,1]) 
            @test_throws ArgumentError diff_action(dpr(-2), [1,-1,1]) 
            
            @test ngens(dpr) == 3
            @test diff_action(f, 1) == dpr[1, [1,0,0]] * u2 + u1 * dpr[2, [1,0,0]]
            @test ngens(dpr) == 5
            @test diff_action(f, 2) == dpr[1, [0,1,0]] * u2 + u1 * dpr[2, [0,1,0]]
            @test ngens(dpr) == 7
            @test diff_action(f, 3) == dpr[1, [0,0,1]] * u2 + u1 * dpr[2, [0,0,1]]
            @test ngens(dpr) == 9

            @test diff_action(g, 1) == -6*dpr[1, [1,0,0]] * u1 * u3 - 3*u1^2 * dpr[3, [1,0,0]] + 4*dpr[2, [1,0,0]]
            @test ngens(dpr) == 10
            @test diff_action(g, 2) == -6*dpr[1, [0,1,0]] * u1 * u3 - 3*u1^2 * dpr[3, [0,1,0]] + 4*dpr[2, [0,1,0]]
            @test ngens(dpr) == 11
            @test diff_action(g, 3) == -6*dpr[1, [0,0,1]] * u1 * u3 - 3*u1^2 * dpr[3, [0,0,1]] + 4*dpr[2, [0,0,1]]
            @test ngens(dpr) == 12

            @test diff_action(f, [2,0,0]) == dpr[1, [2,0,0]] * u2 + 2*dpr[1, [1,0,0]] * dpr[2, [1,0,0]] + u1 * dpr[2, [2,0,0]]
            @test ngens(dpr) == 14
            @test diff_action(f, [0,2,0]) == dpr[1, [0,2,0]] * u2 + 2*dpr[1, [0,1,0]] * dpr[2, [0,1,0]] + u1 * dpr[2, [0,2,0]]
            @test ngens(dpr) == 16
            @test diff_action(f, [0,0,2]) == dpr[1, [0,0,2]] * u2 + 2*dpr[1, [0,0,1]] * dpr[2, [0,0,1]] + u1 * dpr[2, [0,0,2]]
            @test ngens(dpr) == 18

            @test diff_action(g, [1,1,1]) == diff_action(-6*dpr[1, [0,0,1]] * u1 * u3 - 3 * u1^2 * dpr[3, [0,0,1]] + 4*dpr[2, [0,0,1]], [1,1,0])
            @test diff_action(g, [1,1,1]) == diff_action(-6*dpr[1, [0,1,0]] * u1 * u3 - 3 * u1^2 * dpr[3, [0,1,0]] + 4*dpr[2, [0,1,0]], [1,0,1])
            @test diff_action(g, [1,1,1]) == diff_action(-6*dpr[1, [1,0,0]] * u1 * u3 - 3 * u1^2 * dpr[3, [1,0,0]] + 4*dpr[2, [1,0,0]], [0,1,1])
            @test ngens(dpr) == 29  
          end
        end
      end #End for loop
    end #End further constructions
  end #Construction and basic field acces

end #All tests
