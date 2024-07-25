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
          -1 0 0 -1 2]), ["A3 B2", "B2 A3"], :all), # 32
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
