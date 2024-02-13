import Oscar.AbstractAlgebra
import Oscar.AbstractAlgebra: Group

include(
  joinpath(dirname(pathof(AbstractAlgebra)), "..", "test", "Groups-conformance-tests.jl")
)

@testset "LieAlgebras.WeylGroup" begin
  b3_w0 = UInt8[3, 2, 3, 1, 2, 3, 1, 2, 1]
  b4_w0 = UInt8[4, 3, 4, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1]
  f4_w0 = UInt8[4, 3, 2, 3, 1, 2, 3, 4, 3, 2, 3, 1, 2, 3, 4, 3, 2, 3, 1, 2, 3, 1, 2, 1]
  g2_w0 = UInt8[2, 1, 2, 1, 2, 1]

  @testset "weyl_group(::ZZMatrix)" begin
    W = weyl_group(cartan_matrix(:A, 2))
    @test isfinite(W) == true
    @test ngens(W) == 2

    W = weyl_group(cartan_matrix(:A, 3))
    @test isfinite(W) == true
    @test ngens(W) == 3

    W = weyl_group(cartan_matrix(:B, 2))
    @test isfinite(W) == true
    @test ngens(W) == 2

    W = weyl_group(ZZ[2 -2; -2 2]) # TODO: replace with cartan_matrix(A_1^(1)), once functionality for affine type is added
    @test isfinite(W) == false
  end

  @testset "weyl_group(::Symbol, ::Int)" begin
    @test weyl_group(:A, 2) isa WeylGroup
    @test weyl_group(:B, 4) isa WeylGroup
    @test weyl_group(:C, 3) isa WeylGroup
    @test weyl_group(:D, 5) isa WeylGroup
    @test weyl_group(:E, 7) isa WeylGroup
    @test weyl_group(:F, 4) isa WeylGroup
    @test weyl_group(:G, 2) isa WeylGroup

    @test_throws ArgumentError weyl_group(:F, 2)
  end

  @testset "weyl_group(::Tuple{Symbol, Int}...)" begin
    @test weyl_group((:A, 2), (:B, 4)) isa WeylGroup
    @test weyl_group((:C, 3), (:D, 5)) isa WeylGroup
    @test weyl_group((:E, 7)) isa WeylGroup
    @test weyl_group((:F, 4), (:G, 2)) isa WeylGroup

    @test_throws ArgumentError weyl_group((:F, 2), (:B, 4))
    @test_throws ArgumentError weyl_group((:B, 2), (:G, 4))
  end

  @testset "WeylGroup Group conformace test for $(Gname)" for (Gname, G) in [
    ("A1", weyl_group(:A, 1)),
    ("B4", weyl_group(root_system(:B, 4))),
    ("D5", weyl_group(cartan_matrix(:D, 5))),
    ("F4+G2", weyl_group((:F, 4), (:G, 2))),
    ("E6+C3", weyl_group([(:E, 6), (:C, 3)])),
    ("A_1^(1)", weyl_group(ZZ[2 -2; -2 2])), # TODO: replace with cartan_matrix(A_1^(1)), once functionality for affine type is added
    (
      "complicated case 1",
      begin
        cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
        for _ in 1:50
          i, j = rand(1:nrows(cm), 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end
        weyl_group(cm)
      end,
    ),
    (
      "complicated case 2",
      begin
        cm = cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2))
        for _ in 1:50
          i, j = rand(1:nrows(cm), 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end
        weyl_group(root_system(cm))
      end,
    ),
  ]
    # TODO: @felix-roehrich make this work
    # test_Group_interface(G)
    # test_GroupElem_interface(rand(G, 2)...)
  end

  @testset "inv(x::WeylGroupElem)" begin
    W = weyl_group(:A, 2)
    s = gens(W)
    @test inv(s[1]) == s[1]
    @test inv(s[1] * s[2]) == s[2] * s[1]
    @test inv(s[2] * s[1] * s[2]) == s[1] * s[2] * s[1]

    W = weyl_group(:B, 4)
    s = gens(W)
    @test inv(s[2] * s[1]) == s[1] * s[2]
    @test inv(s[3] * s[1]) == s[3] * s[1]
    @test inv(s[2] * s[4] * s[3] * s[4]) == s[4] * s[3] * s[4] * s[2]

    @testset for (fam, rk) in
                 [(:A, 1), (:A, 5), (:B, 3), (:C, 4), (:D, 5), (:F, 4), (:G, 2)]
      W = weyl_group(fam, rk)
      for x in W
        ix = inv(x)
        @test length(ix) == length(x)
        @test isone(ix * x) == isone(x * ix) == true
      end
    end
  end

  @testset "iterate(W::WeylGroup)" begin
    @testset for (fam, rk) in
                 [(:A, 1), (:A, 5), (:B, 3), (:C, 4), (:D, 5), (:F, 4), (:G, 2)]
      W = weyl_group(fam, rk)
      elems = collect(W)
      @test allunique(elems)
      @test length(elems) == order(W)
    end
  end

  @testset "longest_element(W::WeylGroup)" begin
    # A1
    W = weyl_group(:A, 1)
    @test longest_element(W) == gen(W, 1)

    # A2
    W = weyl_group(:A, 2)
    @test word(longest_element(W)) == UInt8[1, 2, 1]

    # B2
    W = weyl_group(:B, 2)
    @test word(longest_element(W)) == UInt8[2, 1, 2, 1]

    # B3
    W = weyl_group(:B, 3)
    @test word(longest_element(W)) == b3_w0

    # F4
    W = weyl_group(:F, 4)
    @test word(longest_element(W)) == f4_w0

    # G2
    W = weyl_group(:G, 2)
    @test word(longest_element(W)) == g2_w0
  end

  @testset "ngens(W::WeylGroup)" begin
    @test ngens(weyl_group(:A, 2)) == 2
    @test ngens(weyl_group(:B, 4)) == 4
    @test ngens(weyl_group(:C, 3)) == 3
    @test ngens(weyl_group(:D, 5)) == 5
    @test ngens(weyl_group(:E, 7)) == 7
    @test ngens(weyl_group(:F, 4)) == 4
    @test ngens(weyl_group(:G, 2)) == 2

    @test ngens(weyl_group((:A, 2), (:B, 4))) == 6
    @test ngens(weyl_group((:C, 3), (:E, 7))) == 10
    @test ngens(weyl_group((:F, 4), (:G, 2))) == 6
  end

  @testset "Base.:(*)(x::WeylGroupElem, y::WeylGroupElem)" begin
    # test short revlex normal form
    W = weyl_group(:A, 2)
    s = gens(W)
    @test parent(s[1] * s[2]) === parent(s[1]) === parent(s[2])

    @test word(s[2] * s[1]) == UInt[2, 1]
    @test word(s[1] * s[2]) == UInt[1, 2]
    @test word(s[1] * s[2] * s[1]) == UInt[1, 2, 1]
    @test word(s[2] * s[1] * s[2]) == UInt[1, 2, 1]

    # test A3
    W = weyl_group(:A, 3)
    s = gens(W)
    @test parent(s[1] * s[2]) === parent(s[1]) === parent(s[2])

    @test word(s[3] * s[1]) == UInt8[3, 1]
    @test word(s[1] * s[3]) == UInt8[3, 1]
    @test word(s[1] * s[3] * s[1]) == UInt8[3]
    @test word(s[3] * s[1] * s[3]) == UInt8[1]
    @test word(s[1] * s[2] * s[1]) == UInt8[1, 2, 1]
    @test word(s[3] * s[2] * s[3]) == UInt8[2, 3, 2]

    # test general multiplication behaviour
    W = weyl_group(:B, 4)
    @test W(b4_w0) == W(b4_w0; normalize=false)

    W = weyl_group(:F, 4)
    @test W(f4_w0) == W(f4_w0; normalize=false)

    W = weyl_group(:G, 2)
    @test W(g2_w0) == W(g2_w0; normalize=false)
  end

  @testset "Base.:(^)(x::WeylGroupElem, n::Int)" begin
    # test A3
    W = weyl_group(:A, 3)
    s = gens(W)

    w = s[1] * s[2] * s[3] * s[2] * s[3]

    @test w^0 == one(W)
    @test w^1 == w
    @test w^2 == w * w
    @test w^3 == w * w * w
    @test w^4 == w * w * w * w
    @test w^5 == w * w * w * w * w
    @test w^6 == w * w * w * w * w * w
    @test w^-1 == inv(w)
    @test w^-2 == inv(w) * inv(w)
    @test w^-3 == inv(w) * inv(w) * inv(w)
    @test w^-4 == inv(w) * inv(w) * inv(w) * inv(w)
  end

  @testset "Base.:(*)(x::WeylGroupElem, w::WeightLatticeElem)" begin
    R = root_system(:A, 2)
    W = weyl_group(R)

    rho = weyl_vector(R)
    @test longest_element(W) * rho == -rho
  end

  @testset "parent(::WeylGroupElem)" begin
    W = weyl_group(:A, 5)
    x = one(W)
    @test parent(x) === x.parent
    @test parent(x) isa WeylGroup

    x = W([1, 3, 5, 4, 2])
    @test parent(x) === x.parent
    @test parent(x) isa WeylGroup
  end

  @testset "ReducedExpressionIterator" begin
    W = weyl_group(:A, 3)
    s = gens(W)

    # test for s1
    iter = reduced_expressions(s[1])
    @test iter.el === s[1]
    @test iter.up_to_commutation == false

    re = collect(iter)
    @test length(re) == 1
    @test re[1] == word(s[1])

    # test for w0
    w0 = longest_element(W)
    iter = reduced_expressions(w0)
    @test iter.el === w0
    @test iter.up_to_commutation == false

    re = collect(iter)
    @test length(re) == 16
    @test re[1] == word(w0)
    @test re[16] == UInt8[3, 2, 1, 3, 2, 3]

    iter = reduced_expressions(w0; up_to_commutation=true)
    @test iter.el === w0
    @test iter.up_to_commutation == true

    re = collect(iter)
    @test length(re) == 8
    @test re[1] == word(w0)
    @test re[8] == UInt8[3, 2, 3, 1, 2, 3]
  end

  @testset "WeylIteratorNoCopy" begin
    WeylIteratorNoCopy = Oscar.LieAlgebras.WeylIteratorNoCopy

    # test simple root systems
    @testset for ((fam, rk), vec) in [
      ((:A, 1), [-42]),
      ((:A, 3), [0, 0, 1]),
      ((:A, 3), [1, 0, 0]),
      ((:A, 5), [1, -1, 2, 0, 2]),
      ((:B, 3), [1, 1, 1]),
      ((:C, 4), [2, 1, 0, 1]),
      ((:D, 5), [-1, 2, 2, -1, -1]),
      ((:E, 6), [1, 2, 0, 0, 2, 1]),
      ((:F, 4), [1, 2, 3, 4]),
      ((:G, 2), [-1, -1]),
    ]
      R = root_system(fam, rk)
      wt = WeightLatticeElem(R, vec)
      dom_wt, conj = conjugate_dominant_weight_with_elem(wt)
      orb = Tuple{WeightLatticeElem,WeylGroupElem}[]
      for tup in WeylIteratorNoCopy(wt)
        push!(orb, deepcopy(tup))
      end

      @test !isnothing(findfirst(==((wt, inv(conj))), orb))
      @test allunique(first.(orb))
      for (ow, x) in orb
        @test x * ow == dom_wt
      end

      gap_num = 0
      gap_W = GAPWrap.WeylGroup(
        GAPWrap.RootSystem(
          GAP.Globals.SimpleLieAlgebra(GAP.Obj(fam), rk, GAP.Globals.Rationals)
        ),
      )
      it = GAPWrap.WeylOrbitIterator(gap_W, GAP.Obj(vec))
      while !GAPWrap.IsDoneIterator(it)
        _ = GAPWrap.NextIterator(it)
        gap_num += 1
      end
      @test length(orb) == gap_num
    end

    # test composite root systems
    @testset for (type, vec) in [
      ([(:A, 1), (:A, 3), (:A, 3)], [-3, 0, 0, 1, 1, 0, 0]),
      ([(:A, 5), (:B, 3)], [1, -1, 2, 0, 2, 1, 1, 1]),
      ([(:C, 2), (:D, 5)], [0, 1, -1, 2, 2, -1, -1]),
      ([(:E, 6)], [1, 2, 0, 0, 2, 1]),
      ([(:F, 4), (:G, 2)], [1, 2, 3, 4, -1, -1]),
    ]
      R = root_system(type...)
      wt = WeightLatticeElem(R, vec)
      dom_wt, conj = conjugate_dominant_weight_with_elem(wt)
      orb = Tuple{WeightLatticeElem,WeylGroupElem}[]
      for tup in WeylIteratorNoCopy(wt)
        push!(orb, deepcopy(tup))
      end

      @test !isnothing(findfirst(==((wt, inv(conj))), orb))
      @test allunique(first.(orb))
      for (ow, x) in orb
        @test x * ow == dom_wt
      end

      gap_num = 0
      gap_L = GAP.Globals.DirectSumOfAlgebras(
        GAP.Obj([
          GAP.Globals.SimpleLieAlgebra(GAP.Obj(fam), rk, GAP.Globals.Rationals) for
          (fam, rk) in type
        ]),
      )
      gap_W = GAPWrap.WeylGroup(GAPWrap.RootSystem(gap_L))
      it = GAPWrap.WeylOrbitIterator(gap_W, GAP.Obj(vec))
      while !GAPWrap.IsDoneIterator(it)
        _ = GAPWrap.NextIterator(it)
        gap_num += 1
      end
      @test length(orb) == gap_num
    end
  end

  @testset "WeylOrbitIterator" begin
    @test eltype(WeylOrbitIterator) == WeightLatticeElem

    @testset for ((fam, rk), vec) in [
      ((:A, 1), [-42]),
      ((:A, 3), [0, 0, 1]),
      ((:A, 3), [1, 0, 0]),
      ((:A, 5), [1, -1, 2, 0, 2]),
      ((:B, 3), [1, 1, 1]),
      ((:C, 4), [2, 1, 0, 1]),
      ((:D, 5), [-1, 2, 2, -1, -1]),
      ((:E, 6), [1, 2, 0, 0, 2, 1]),
      ((:F, 4), [1, 2, 3, 4]),
      ((:G, 2), [-1, -1]),
    ]
      R = root_system(fam, rk)
      wt = WeightLatticeElem(R, vec)
      dom_wt, conj = conjugate_dominant_weight_with_elem(wt)
      orb = collect(WeylOrbitIterator(wt))

      @test !isnothing(findfirst(==(wt), orb))
      @test allunique(orb)
    end
  end
end
