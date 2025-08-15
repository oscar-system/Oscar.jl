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
  __update_internals! = Oscar.__update_internals!
  
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

        @testset "Check public fields at construction" begin
          @test base_ring(dpr) == ZZ
          @test all(var -> base_ring(var) == ZZ, vars)
          @test elementary_symbols(dpr) == [:u1, :u2, :u3]
          @test ndiffs(dpr) == 3
          @test nelementary_symbols(dpr) == 3
          @test all(var -> parent(var) === dpr, vars)

          ran = ranking(dpr)
          if dpr isa DifferencePolyRing
            @test ran isa DifferenceRanking
          end
          if dpr isa DifferentialPolyRing
            @test ran isa DifferentialRanking
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
            @test nelementary_symbols(dpr) == 3

            ran = ranking(dpr)
            if dpr isa DifferencePolyRing
              @test ran isa DifferenceRanking
            end
            if dpr isa DifferentialPolyRing
              @test ran isa DifferentialRanking
            end
            @test partition(ran) == [[1,0,0], [0,1,0], [0,0,1]]
            @test index_ordering_matrix(ran) == ZZ[0 0 1; 0 1 0; 1 0 0]
            @test riquier_matrix(ran) == ZZ[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0; 0 0 0 1 0 0]
          end

          @testset "Some computations" begin
            @test dpr() == zero(dpr)
            @test dpr(0) == zero(dpr)
            @test is_zero(dpr())
            @test !is_unit(dpr())
            @test is_constant(dpr())
            @test_throws ErrorException inv(dpr())
            @test trailing_coefficient(dpr()) == zero(ZZ)
            @test leading_coefficient(dpr()) == zero(ZZ)
            @test_throws ArgumentError leading_term(dpr(0))
            @test_throws ArgumentError leading_monomial(dpr(0))
            @test length(dpr()) == 0

            @test dpr(1) == one(dpr)
            @test is_one(dpr(1))
            @test is_unit(dpr(1))
            @test is_constant(dpr())
            @test inv(dpr(1)) == 1
            @test trailing_coefficient(dpr(1)) == one(ZZ)
            @test leading_coefficient(dpr(1)) == one(ZZ)
            @test leading_term(dpr(1)) == dpr(1)
            @test leading_monomial(dpr(1)) == dpr(1)
            @test length(dpr(1)) == 1

            @test !is_unit(dpr(-2))
            @test is_constant(dpr(-2))
            @test_throws ErrorException inv(dpr(-2))
            @test trailing_coefficient(dpr(-2)) == ZZ(-2)
            @test leading_coefficient(dpr(-2)) == ZZ(-2)
            @test leading_term(dpr(-2)) == dpr(-2)
            @test leading_monomial(dpr(-2)) == dpr(1)
            @test length(dpr(-2)) == 1
          end
        
        end #End check adding variables

      end #End for loop

    end #Construction with keywords

  end #Construction and basic field acces

end #All tests
