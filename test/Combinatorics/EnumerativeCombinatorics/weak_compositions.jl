@testset "weak compositions for integer type $T" for T in [Int, Int8, ZZRingElem]
  @test length(weak_compositions(T(3), T(2))) == 4
  @test (@inferred collect(weak_compositions(T(3), T(2)))) == map(weak_composition, [T[3, 0], T[2, 1], T[1, 2], T[0, 3]])
  @test isempty(weak_compositions(T(1), T(0)))
  @test (@inferred collect(weak_compositions(T(1), T(0)))) == Oscar.WeakComposition{T}[]
  @test length(weak_compositions(T(0), T(0))) == 1
  @test (@inferred collect(weak_compositions(T(0), T(0)))) == [weak_composition(T[])]
  @test length(weak_compositions(T(0), T(3))) == 1
  @test (@inferred collect(weak_compositions(T(0), T(3)))) == [weak_composition(T[0, 0, 0])]

  @testset "inplace iteration" begin
    Ci = weak_compositions(T(5),T(3), inplace=true)
    Cf = weak_compositions(T(5),T(3))
    @test all(splat(==), zip(Ci, Cf))
  end

  l = @inferred collect(weak_compositions(T(5), T(3)))
  @test length(l) == 21
  @test all(x -> sum(x) == 5, l)
  b = unique!(l)
  @test length(b) == length(l)

  l = @inferred collect(weak_compositions(T(5), T(0)))
  @test l == Vector{Oscar.WeakComposition{T}}()
  l = @inferred collect(weak_compositions(T(0), T(7)))
  @test l == [weak_composition(T[0 for i in 1:7])]
  l = @inferred collect(weak_compositions(T(1), T(7)))
  @test length(l) == 7
  @test allunique(l)

  l = @inferred collect(weak_compositions(T(7), T(1)))
  @test l == [weak_composition(T[7])]

  @test number_of_weak_compositions(T(1), T(-2)) == ZZ(0)
  @test number_of_weak_compositions(T(-1), T(0)) == ZZ(0)
end
