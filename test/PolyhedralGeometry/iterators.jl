@testset "iterators" begin
  @testset "SubObjectIterator" begin
    # SubObjectIterator is flexible due to the `Function` object it stores, which is used for `getindex`.
    # For the tests we want this iterator to return a `Pair` consisting of a vertex and a halfspace.
    # This context is very specific and probably not meaningful from a mathematical point of view,
    # but it serves well for abstractly testing this type.
    c = cube(2)
    function _testSOI(
      ::Type{Pair{Matrix{QQFieldElem},QQFieldElem}},
      obj::Polyhedron{QQFieldElem},
      i::Base.Integer,
    )
      x = Oscar.pm_object(obj).VERTICES[i, 2:end]
      # x = PointVector{QQFieldElem}(obj.VERTICES[i, 2:end])
      a, b = Oscar.decompose_hdata(Oscar.pm_object(obj).FACETS)
      # y = AffineHalfspace{QQFieldElem}(a[i, :], b[i])
      return Pair{Matrix{QQFieldElem},QQFieldElem}(
        convert(Matrix{QQFieldElem}, hcat(x, x)), b[i]
      )
    end
    soi = SubObjectIterator{Pair{Matrix{QQFieldElem},QQFieldElem}}(c, _testSOI, 4)
    @test soi isa AbstractVector{Pair{Matrix{QQFieldElem},QQFieldElem}}
    @test length(soi) == 4
    @test firstindex(soi) == 1
    @test lastindex(soi) == 4
    for i in 1:4
      p = soi[i]
      @test p[1] == convert(
        Matrix{QQFieldElem},
        hcat(Oscar.pm_object(c).VERTICES[i, 2:end], Oscar.pm_object(c).VERTICES[i, 2:end]),
      )
      a, b = Oscar.decompose_hdata(Oscar.pm_object(c).FACETS)
      @test p[2] == convert(QQFieldElem, b[i])
    end
    @test_throws ArgumentError ray_indices(soi)
    @test_throws ArgumentError linear_inequality_matrix(soi)
    @test_throws ArgumentError halfspace_matrix_pair(soi)
    @test_throws ArgumentError Oscar.linear_matrix_for_polymake(soi)
    @test_throws ArgumentError Oscar.affine_matrix_for_polymake(soi)
    @test_throws ArgumentError Oscar.matrix_for_polymake(soi)
    soi2 = SubObjectIterator{PointVector{QQFieldElem}}(c, _testSOI, 4)
    @test_throws ArgumentError point_matrix(soi2)
  end
  @testset "RayVector" begin
    rv = ray_vector([1, 0, 0])
    F = facets(cube(3))
    @test_throws ArgumentError rv in F[1]
  end
end
