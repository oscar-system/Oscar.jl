@testset "LieAlgebras.AbstractLieAlgebra" begin
  function sl2_struct_consts(R::Field)
    sc = zeros(R, 3, 3, 3)
    sc[1, 2, 3] = R(1)
    sc[2, 1, 3] = R(-1)
    sc[3, 1, 1] = R(2)
    sc[1, 3, 1] = R(-2)
    sc[3, 2, 2] = R(-2)
    sc[2, 3, 2] = R(2)
    return sc
  end

  function G2_non_canonical_basis_struct_consts(R::Field)
    # shared by Willem de Graaf on slack
    tab = fill(sparse_row(R), 14, 14)
    #! format: off
    tab[1,3] = sparse_row(R, [3,8,10,11,13], [3,2,-4,-2,2])
    tab[1,4] = sparse_row(R, [4,8,9,10,11,13], [5,4//3,-3,-8//3,8//3,-8//3])
    tab[1,5] = sparse_row(R, [4,5,8,9], [-1,-2,1,1])
    tab[1,6] = sparse_row(R, [4,5,6,8,9,10,11,13], [-1,-5,3,3,1,-1,1,-1])
    tab[1,7] = sparse_row(R, [7,8,10,11,13], [-1,-2//3,4//3,-4//3,4//3])
    tab[1,8] = sparse_row(R, [8,10,11,13], [-1//3,2//3,-2//3,2//3])
    tab[1,9] = sparse_row(R, [4,8,9,10,11,13], [6,2,-4,-4,4,-4])
    tab[1,10] = sparse_row(R, [8,10,11,13], [2//3,-4//3,-5//3,5//3])
    tab[1,11] = sparse_row(R, [4,8,9,10,11,13], [1,4//3,-1,-8//3,-1//3,1//3])
    tab[1,12] = sparse_row(R, [8,10,11,12,13], [2,-1,1,-3,-1])
    tab[1,13] = sparse_row(R, [4,9], [1,-1])
    tab[1,14] = sparse_row(R, [8,10,11,12,13,14], [2//3,5//3,7//3,-4,-7//3,1])
    tab[2,3] = sparse_row(R, [3,8,10,11,13], [-2,-4//3,8//3,4//3,-4//3])
    tab[2,4] = sparse_row(R, [4,8,9,10,11,13], [-2,-1//3,1,2//3,-2//3,2//3])
    tab[2,5] = sparse_row(R, [4,5,8,9,10,11,13], [1,1,1//3,-1,-2//3,2//3,-2//3])
    tab[2,6] = sparse_row(R, [4,5,6,8,9], [1,2,-1,-1,-1])
    tab[2,7] = sparse_row(R, [7,8,10,11,13], [1,2//3,-4//3,4//3,-4//3])
    tab[2,8] = sparse_row(R, [8,10,11,13], [5//3,-4//3,4//3,-4//3])
    tab[2,9] = sparse_row(R, [4,8,9,10,11,13], [-2,-2//3,1,4//3,-4//3,4//3])
    tab[2,10] = sparse_row(R, [8,10,11,13], [1//3,1//3,5//3,-5//3])
    tab[2,11] = sparse_row(R, [4,8,9,10,11,13], [-1,-4//3,1,8//3,-2//3,-1//3])
    tab[2,12] = sparse_row(R, [8,10,11,12,13], [-1//3,2//3,1//3,1,-1//3])
    tab[2,13] = sparse_row(R, [4,8,9,10,11,13], [-1,-1//3,1,2//3,-2//3,-1//3])
    tab[2,14] = sparse_row(R, [8,10,11,12,13], [1//3,-2//3,-1//3,1,1//3])
    tab[3,4] = sparse_row(R, [4,5,7,8,9], [-1,-1,1//2,1//2,1])
    tab[3,5] = sparse_row(R, [5,6,8,10,11,13], [2,-2,-2//3,-2//3,2//3,-2//3])
    tab[3,6] = sparse_row(R, [10,11,13], [-1,1,-1])
    tab[3,7] = sparse_row(R, [4,5,8,9,10,11,13], [7,-1,3,-3,-5,5,-5])
    tab[3,8] = sparse_row(R, [4,5,6,8,9,10,11,13], [-1,3,-4,-7//3,1,5//3,-5//3,5//3])
    tab[3,9] = sparse_row(R, [7,8,10,11,13], [1//2,1//6,-1//3,1//3,-1//3])
    tab[3,10] = sparse_row(R, [2,4,5,6,8,9,10,11,13], [-2,-1,1,-2,-1,1,1,-1,1])
    tab[3,11] = sparse_row(R, [2,8,10,11,12,13], [-2,-1//3,-1//3,-2//3,1,2//3])
    tab[3,12] = sparse_row(R, [2,4,5,6,9,10,11], [-2,2,2,-2,-2,-2,2])
    tab[3,14] = sparse_row(R, [4,5,6,9,10,11], [2,2,-2,-2,-2,2])
    tab[4,5] = sparse_row(R, [1,4,7,8,9,10,11,12,13,14], [1,-1,-1//2,-3//2,1,6,-2,-4,-1,1])
    tab[3,13] = sparse_row(R, [4,5,9,11,12,13], [-1,-1,1,-1,1,1])
    tab[4,6] = sparse_row(R, [1,4,8,10,11,12,13], [1,1,-1//3,11//3,-2//3,-3,-7//3])
    tab[4,7] = sparse_row(R, [1,2,3,4,5,8,9,10,11,12,13,14], [-2,-6,3//2,7,4,5//2,-7,-9,7//2,4,-1//2,-4])
    tab[4,8] = sparse_row(R, [3,4,7,8,9,10,11,12,13,14], [3//2,3,-1,3//2,-3,-3,7//2,-2,-1//2,2])
    tab[4,9] = sparse_row(R, [3,4,8,9,10,11,13], [-3//2,-6,-31//6,6,31//3,-53//6,17//6])
    tab[4,10] = sparse_row(R, [3,4,5,7,8,9,10,11,12,13,14], [3//2,4,1,-1,3//2,-4,-4,7//2,-1,-1//2,1])
    tab[4,11] = sparse_row(R, [4,5,7,8,9,10,11,13], [-2,1,-1//2,-19//6,2,16//3,-16//3,7//3])
    tab[4,12] = sparse_row(R, [5,7,8,10,11,12,13,14], [1,-1,-4//3,5//3,-2//3,-1,2//3,1])
    tab[4,13] = sparse_row(R, [3,8,10,11,13], [3//2,-1//6,1//3,-11//6,11//6])
    tab[4,14] = sparse_row(R, [1,2,4,5,6,7,8,9,11,12,13,14], [-2,-3,3,3,-3,-1//2,-3//2,-1,1,-1,-1,1])
    tab[5,6] = sparse_row(R, [4,7,8,9,10,11,12,13,14], [2,1//2,7//6,-1,-7//3,4//3,1,-4//3,-1])
    tab[5,7] = sparse_row(R, [4,5,8,9,10,11,12,13,14], [5,4,7//3,-5,-26//3,5//3,1,4//3,-1])
    tab[5,8] = sparse_row(R, [4,7,8,9,10,11,12,13,14], [1,-1,4//3,-1,-8//3,5//3,1,4//3,-1])
    tab[5,9] = sparse_row(R, [1,7,8,10,11,12,13], [-1,1,-1,-4,-2,6,2])
    tab[5,10] = sparse_row(R, [4,7,8,9,10,11,12,13,14], [1,-1//2,3//2,-1,-3,2,1,1,-1])
    tab[5,11] = sparse_row(R, [1,2,7,8,10,11,12,13], [1,2,1//2,-5//6,-4//3,-5//3,3,5//3])
    tab[5,12] = sparse_row(R, [7,8,10,11,13], [-1//2,1//6,-1//3,-2//3,2//3])
    tab[5,13] = sparse_row(R, [1,2,4,7,8,9,10,11,12,13,14], [1,2,1,1//2,5//6,-1,-14//3,2//3,4,7//3,-1])
    tab[5,14] = sparse_row(R, [1,2,7,8,10,11,13], [-2,-3,-3//2,-1//6,1//3,-4//3,4//3])
    tab[6,7] = sparse_row(R, [4,5,8,9,10,11,13], [5,4,2,-5,-8,2,1])
    tab[6,8] = sparse_row(R, [4,7,8,9,10,11,13], [1,-1,1,-1,-2,2,1])
    tab[6,9] = sparse_row(R, [1,4,8,9,10,11,12,13], [-1,-4,-8//3,2,-2//3,-16//3,6,16//3])
    tab[6,10] = sparse_row(R, [4,7,8,9,10,11,13], [1,-1//2,11//6,-1,-8//3,8//3,1//3])
    tab[6,11] = sparse_row(R, [3,4,8,9,10,11,12,13], [-1//2,-2,-7//6,1,1//3,-17//6,3,17//6])
    tab[6,12] = sparse_row(R, [1,2,7,8], [1,1,-1//2,1//2])
    tab[6,13] = sparse_row(R, [3,4,8,10,11,12,13], [-1//2,-1,-1//2,-2,-1//2,3,7//2])
    tab[6,14] = sparse_row(R, [1,2,7,8,10,11,13], [-1,-2,-3//2,-1//2,1,-1,1])
    tab[7,8] = sparse_row(R, [1,2,8,10,11,12,13,14], [2,6,-2//3,4//3,2//3,-2,-2//3,2])
    tab[7,9] = sparse_row(R, [3,4,5,8,9,10,11,12,13,14], [-3//2,-14,-8,-19//6,14,43//3,-53//6,-4,17//6,4])
    tab[7,10] = sparse_row(R, [1,2,4,5,8,9,10,11,12,13,14], [2,6,1,1,-2//3,-1,1//3,2//3,-1,-2//3,1])
    tab[7,11] = sparse_row(R, [4,5,8,9,10,11,13], [-4,-3,-1//3,4,11//3,-11//3,2//3])
    tab[7,12] = sparse_row(R, [4,5,8,9,10,11,12,13,14], [1,1,-2//3,-1,1//3,2//3,-1,-2//3,1])
    tab[7,13] = sparse_row(R, [1,2,4,5,9,10,11,13], [2,6,-5,-4,5,4,-4,1])
    tab[7,14] = sparse_row(R, [4,8,9,10,11,12,13,14], [4,5,-2,-4,5,-1,-5,1])
    tab[8,9] = sparse_row(R, [3,4,7,8,9,10,11,13], [-3//2,-6,2,-23//6,6,23//3,-37//6,1//6])
    tab[8,10] = sparse_row(R, [4,5,9,10,12,14], [1,1,-1,-1,1,-1])
    tab[8,11] = sparse_row(R, [1,2,4,5,7,8,9,10,11], [2,4,-2,1,1,-2,2,3,-3])
    tab[8,12] = sparse_row(R, [4,5,9,10,12,14], [1,1,-1,-1,-1,1])
    tab[8,13] = sparse_row(R, [1,2,4,7,8,9,10,11,12,13,14], [2,4,-3,1,-1,3,2,-4,2,1,-2])
    tab[8,14] = sparse_row(R, [4,8,9,10,11,12,13,14], [4,5//3,-2,-10//3,7//3,-1,-7//3,1])
    tab[9,10] = sparse_row(R, [3,4,7,8,9,10,11,13], [3//2,6,-3//2,4,-6,-8,13//2,-1//2])
    tab[9,11] = sparse_row(R, [7,8,10,11,13], [-1//2,-5//6,5//3,-5//3,5//3])
    tab[9,12] = sparse_row(R, [4,7,8,9,10,11,13], [-1,-3//2,-5//6,1,5//3,-5//3,5//3])
    tab[9,13] = sparse_row(R, [3,4,8,9,10,11,13], [3//2,6,23//6,-6,-23//3,37//6,-1//6])
    tab[9,14] = sparse_row(R, [1,2,4,5,6,7,8,9,10,11,13], [-4,-6,-1,3,-3,-1,-8//3,1,7//3,-7//3,7//3])
    tab[10,11] = sparse_row(R, [1,2,4,5,7,8,9,10,11,12,13], [1,2,-2,1,1//2,-5//2,2,3,-4,1,1])
    tab[10,12] = sparse_row(R, [4,5,8,9,10,11,12,13,14], [1,1,-1//3,-1,-1//3,1//3,-1,-1//3,1])
    tab[10,13] = sparse_row(R, [1,2,4,5,7,8,9,10,11,12,13,14], [1,2,-4,-1,1//2,-3//2,4,3,-5,2,2,-1])
    tab[10,14] = sparse_row(R, [4,8,9,10,11,12,13,14], [4,4//3,-2,-8//3,8//3,-1,-8//3,1])
    tab[11,12] = sparse_row(R, [1,2,7,8,10,11,12,13], [-1,-2,-1//2,1//6,2//3,1//3,-1,-1//3])
    tab[11,13] = sparse_row(R, [4,5,8,9,10,11,12], [2,-1,2,-2,-4,3,1])
    tab[11,14] = sparse_row(R, [1,2,7,8,10,11,13], [-3,-5,-1//2,1//6,-1//3,1//3,-1//3])
    tab[12,13] = sparse_row(R, [1,2,4,5,7,8,9,10,11,12,13,14], [1,2,-1,-1,1//2,1//2,1,-1,-1,2,1,-1])
    tab[12,14] = sparse_row(R, [4,5,8,9,10,11,13], [1,1,-1//3,-1,-1//3,1//3,-1//3])
    tab[13,14] = sparse_row(R, [1,2,4,7,8,9,10,11,12,13,14], [-3,-5,4,-1//2,7//6,-2,-7//3,10//3,-1,-10//3,1])
    #! format: on
    tab -= permutedims(tab)
    return tab
  end

  @testset "conformance tests" begin
    @testset "0-dim Lie algebra /QQ" begin
      L = lie_algebra(QQ, Matrix{sparse_row_type(QQFieldElem)}(undef, 0, 0), Symbol[])
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
      @test is_abelian(L)
    end

    @testset "sl_2(QQ) using structure constants" begin
      L = lie_algebra(QQ, sl2_struct_consts(QQ), ["e", "f", "h"])
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "sl_2(CF(4)) using structure constants" begin
      CF4 = cyclotomic_field(4)[1]
      L = lie_algebra(CF4, sl2_struct_consts(CF4), ["e", "f", "h"])
      lie_algebra_conformance_test(
        L,
        AbstractLieAlgebra{AbsSimpleNumFieldElem},
        AbstractLieAlgebraElem{AbsSimpleNumFieldElem},
      )
    end

    @testset "4-dim abelian Lie algebra /QQ" begin
      L = abelian_lie_algebra(AbstractLieAlgebra, QQ, 4)
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
      @test is_abelian(L)
    end

    @testset "A_4(QQ)" begin
      L = lie_algebra(QQ, :A, 4)
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "B_3(QQ)" begin
      L = lie_algebra(QQ, :B, 3)
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "A_4(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], :A, 4)
      lie_algebra_conformance_test(
        L,
        AbstractLieAlgebra{AbsSimpleNumFieldElem},
        AbstractLieAlgebraElem{AbsSimpleNumFieldElem},
      )
    end

    @testset "B_3(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], :B, 3)
      lie_algebra_conformance_test(
        L,
        AbstractLieAlgebra{AbsSimpleNumFieldElem},
        AbstractLieAlgebraElem{AbsSimpleNumFieldElem},
      )
    end

    @testset "G_2(QQ) with non-canonical basis" begin
      L = lie_algebra(QQ, G2_non_canonical_basis_struct_consts(QQ))
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "G_2(GF(137)) with non-canonical basis" begin
      F = GF(137)
      L = lie_algebra(F, G2_non_canonical_basis_struct_consts(F))
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{FqFieldElem}, AbstractLieAlgebraElem{FqFieldElem}
      )
    end
  end

  @testset "Lie algebra from root system" begin
    function test_lie_algebras_from_root_system(
      rs::RootSystem;
      repeats::Union{Oscar.IntegerUnion,Symbol}=:all,
      gap_type::Union{Vector{String},Nothing}=nothing,
    )
      test_lie_algebras_from_root_system(QQ, rs; repeats, gap_type)
    end

    function test_lie_algebras_from_root_system(
      R::Field,
      rs::RootSystem;
      repeats::Union{Oscar.IntegerUnion,Symbol}=:all,
      gap_type::Union{Vector{String},Nothing}=nothing,
    )
      @req (repeats isa Oscar.IntegerUnion && repeats >= 0) || repeats == :all "`repeats` must be a positive integer or :all"
      n_freedom_degs = n_positive_roots(rs) - n_simple_roots(rs)
      if n_freedom_degs == 0
        signs = Bool[]
        test_lie_algebra_from_root_system(R, rs, signs; gap_type)
      elseif repeats == :all
        for signs in AbstractAlgebra.ProductIterator([true, false], n_freedom_degs)
          test_lie_algebra_from_root_system(R, rs, signs; gap_type)
        end
      else
        for _ in 1:repeats
          signs = rand(Bool, n_freedom_degs)
          test_lie_algebra_from_root_system(R, rs, signs; gap_type)
        end
      end
      return nothing
    end

    function test_lie_algebra_from_root_system(
      R::Field,
      rs::RootSystem,
      signs::Vector{Bool};
      gap_type::Union{Vector{String},Nothing}=nothing,
    )
      struct_consts = Oscar.LieAlgebras._struct_consts(R, rs, signs)
      # `check=true`` is slow but tests if the structure constants define a Lie algebra, in particular if they satisfy the Jacobi identity
      L = lie_algebra(R, struct_consts; check=true)
      @test true # count the number of testcases
      if !isnothing(gap_type)
        @test String(
          GAP.Globals.SemiSimpleType(
            codomain(Oscar._iso_oscar_gap(L; set_attributes=false)) # let GAP infer the type
          ),
        ) in gap_type
      end
    end

    testcases_fast = Tuple{
      String,RootSystem,Vector{String},Union{Oscar.IntegerUnion,Symbol}
    }[
      ("A1", root_system(:A, 1), ["A1"], :all), # 1
      ("A2", root_system(:A, 2), ["A2"], :all), # 2
      ("A3", root_system(:A, 3), ["A3"], :all), # 8
      ("A4", root_system(:A, 4), ["A4"], 5),
      ("A5", root_system(:A, 5), ["A5"], 1),
      ("B2", root_system(:B, 2), ["B2", "C2"], :all), # 4
      ("B3", root_system(:B, 3), ["B3"], 5),
      ("B4", root_system(:B, 4), ["B4"], 1),
      ("C2", root_system(:C, 2), ["B2", "C2"], :all), # 4
      ("C3", root_system(:C, 3), ["C3"], 5),
      ("C4", root_system(:C, 4), ["C4"], 1),
      ("D4", root_system(:D, 4), ["D4"], 8),
      ("D5", root_system(:D, 5), ["D5"], 1),
      ("F4", root_system(:F, 4), ["F4"], 1),
      ("G2", root_system(:G, 2), ["G2"], :all), # 16
      ("A2 + B2", root_system((:A, 2), (:B, 2)), ["A2 B2", "B2 A2"], :all), # 8
      ("A3 + B2, shuffled",
        root_system([
          2 0 0 0 -1;
          0 2 -2 0 0;
          0 -1 2 0 0;
          0 0 0 2 -1;
          -1 0 0 -1 2]), ["A3 B2", "B2 A3"], 5),
    ]

    testcases_indepth = Tuple{
      String,RootSystem,Vector{String},Union{Oscar.IntegerUnion,Symbol}
    }[
      ("A1", root_system(:A, 1), ["A1"], :all), # 1
      ("A2", root_system(:A, 2), ["A2"], :all), # 2
      ("A3", root_system(:A, 3), ["A3"], :all), # 8
      ("A4", root_system(:A, 4), ["A4"], :all), # 64
      ("A5", root_system(:A, 5), ["A5"], 10),
      ("A6", root_system(:A, 6), ["A6"], 5),
      ("A7", root_system(:A, 7), ["A7"], 3),
      ("B2", root_system(:B, 2), ["B2", "C2"], :all), # 4
      ("B3", root_system(:B, 3), ["B3"], :all), # 64
      ("B4", root_system(:B, 4), ["B4"], 20),
      ("B5", root_system(:B, 5), ["B5"], 5),
      ("C2", root_system(:C, 2), ["B2", "C2"], :all), # 4
      ("C3", root_system(:C, 3), ["C3"], :all), # 64
      ("C4", root_system(:C, 4), ["C4"], 20),
      ("C5", root_system(:C, 5), ["C5"], 5),
      ("D4", root_system(:D, 4), ["D4"], 50),
      ("D5", root_system(:D, 5), ["D5"], 5),
      ("D6", root_system(:D, 6), ["D6"], 2),
      ("E6", root_system(:E, 6), ["E6"], 1),
      ("F4", root_system(:F, 4), ["F4"], 5),
      ("F4 again", root_system(transpose(cartan_matrix(:F, 4))), ["F4"], 3),
      ("G2", root_system(:G, 2), ["G2"], :all), # 16
      ("A2 + B3", root_system((:A, 2), (:B, 3)), ["A2 B3", "B3 A2"], 50),
      ("A3 + B2, shuffled",
        root_system([
          2 0 0 0 -1;
          0 2 -2 0 0;
          0 -1 2 0 0;
          0 0 0 2 -1;
          -1 0 0 -1 2]), ["A3 B2", "B2 A3", "A3 C2", "C2 A3"], :all), # 32
      ("B3 + G2, shuffled",
        root_system([
          2 -2 0 0 0;
          -1 2 0 -1 0;
          0 0 2 0 -1;
          0 -1 0 2 0;
          0 0 -3 0 2
        ]), ["B3 G2", "G2 B3"], 20),
    ]

    # After changing the function, verify that the tests still pass with `testcases_indepth`,
    # but this is too slow for the CI
    for (name, rs, gap_type, repeats) in testcases_fast
      @testset "$name" begin
        test_lie_algebras_from_root_system(rs; repeats, gap_type)
      end
    end

    @testset "other fields" begin
      @testset "$R" for R in
                        [cyclotomic_field(4)[1], GF(3), GF(2, 3), algebraic_closure(QQ)]
        test_lie_algebras_from_root_system(R, root_system(:A, 3); repeats=1)
        test_lie_algebras_from_root_system(R, root_system(:B, 3); repeats=1)
        test_lie_algebras_from_root_system(R, root_system(:C, 3); repeats=1)
        test_lie_algebras_from_root_system(R, root_system(:D, 4); repeats=1)
        test_lie_algebras_from_root_system(R, root_system(:F, 4); repeats=1)
        test_lie_algebras_from_root_system(R, root_system(:G, 2); repeats=1)
      end
    end
  end
end
