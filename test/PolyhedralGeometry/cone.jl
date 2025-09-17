@testset "Cone{$T}" for (f, T) in _prepare_scalar_types()
  pts = [1 0; 0 0; 0 1]
  Cone1 = positive_hull(f, pts)
  R = [1 0 0; 0 0 1]
  L = [0 1 0]
  Cone2 = positive_hull(f, R, L)
  Cone3 = positive_hull(f, R, L; non_redundant=true)
  Cone4 = positive_hull(f, R)
  Cone5 = positive_hull(f, [1 0 0; 1 1 0; 1 1 1; 1 0 1])
  Cone6 = positive_hull(f, [1//3 1//2; 4//5 2])
  Cone7 = positive_hull(f, [0 1])
  Cone8 = positive_hull(f, [1 1; 1 -1])

  @testset "core functionality" begin
    @test is_pointed(Cone1)
    @test issubset(Cone7, Cone1)
    @test !issubset(Cone1, Cone7)
    @test [1, 0] in Cone1
    @test !([-1, -1] in Cone1)
    if T == QQFieldElem
      @test !is_smooth(Cone2)
      @test is_smooth(Cone7)
      @test !is_smooth(Cone8)
    end
    @test is_simplicial(Cone7)
    @test !is_simplicial(Cone5)
    @test is_fulldimensional(Cone1)
    if T == QQFieldElem
      @test hilbert_basis(Cone1) isa SubObjectIterator{PointVector{ZZRingElem}}
      @test length(hilbert_basis(Cone1)) == 2
      @test issetequal(hilbert_basis(Cone1), point_vector.(Ref(ZZ), [[1, 0], [0, 1]]))
      @test generator_matrix(hilbert_basis(Cone1)) ==
        _oscar_matrix_from_property(ZZ, hilbert_basis(Cone1))
    end
    @test n_rays(Cone1) == 2
    @test rays(RayVector{T}, Cone1) isa SubObjectIterator{RayVector{T}}
    @test rays(Cone1) isa SubObjectIterator{RayVector{T}}
    @test rays(RayVector, Cone1) isa SubObjectIterator{RayVector{T}}
    @test issetequal(rays(Cone1), ray_vector.(Ref(f), [[1, 0], [0, 1]]))
    @test vector_matrix(rays(Cone1)) == _oscar_matrix_from_property(f, rays(Cone1))
    if T == QQFieldElem
      @test matrix(QQ, rays(Cone1)) == _oscar_matrix_from_property(f, rays(Cone1))
      let r = rays(Cone6)
        m = matrix(ZZ, r[1] == [1//3, 1//2] ? [2 3; 2 5] : [2 5; 2 3])
        @test matrix(ZZ, rays(Cone6)) == m
      end
    end
    @test length(rays(Cone1)) == 2
    for S in [LinearHalfspace{T}, Cone{T}]
      @test facets(S, Cone1) isa SubObjectIterator{S}
      @test length(facets(S, Cone1)) == 2
      if S == LinearHalfspace{T}
        @test issetequal(facets(S, Cone1), linear_halfspace.(Ref(f), [[-1, 0], [0, -1]]))
      else
        @test issetequal(facets(S, Cone1), positive_hull.(Ref(f), [[1 0], [0 1]]))
      end
      @test linear_inequality_matrix(facets(S, Cone1)) ==
        _oscar_matrix_from_property(f, facets(S, Cone1))
      @test Oscar.linear_matrix_for_polymake(facets(S, Cone1)) ==
        _polymake_matrix_from_property(facets(S, Cone1))
      @test _check_im_perm_rows(ray_indices(facets(S, Cone1)), [[1], [2]])
      @test _check_im_perm_rows(incidence_matrix(facets(S, Cone1)), [[1], [2]])
    end
    @test _check_im_perm_rows(facets(IncidenceMatrix, Cone1), [[1], [2]])
    @test facets(Halfspace, Cone1) isa SubObjectIterator{LinearHalfspace{T}}
    @test facets(Cone1) isa SubObjectIterator{LinearHalfspace{T}}
    @test linear_span(Cone4) isa SubObjectIterator{LinearHyperplane{T}}
    @test length(linear_span(Cone4)) == 1
    @test linear_span(Cone4) == [linear_hyperplane(f, [0 1 0])]
    @test linear_equation_matrix(linear_span(Cone4)) == matrix(f, [0 1 0])

    @test !is_pointed(Cone2)
    @test !is_pointed(Cone3)
    @test !is_fulldimensional(Cone4)
    @test is_fulldimensional(Cone2)
    @test Cone2 == Cone3
    @test Cone4 != Cone2

    @test length(unique([Cone2, Cone3, Cone4])) == 2

    @test dim(Cone4) == 2
    @test dim(Cone2) == 3
    @test ambient_dim(Cone2) == 3
    @test lineality_space(Cone2) isa SubObjectIterator{RayVector{T}}
    @test generator_matrix(lineality_space(Cone2)) == matrix(f, L)
    if T == QQFieldElem
      @test matrix(QQ, lineality_space(Cone2)) == matrix(QQ, L)
      @test matrix(ZZ, lineality_space(Cone2)) == matrix(ZZ, L)
    end
    @test length(lineality_space(Cone2)) == 1
    @test lineality_space(Cone2) == [L[1, :]]
    @test vector_matrix(rays(Cone4)) == _oscar_matrix_from_property(f, rays(Cone4))
    @test codim(Cone4) == 1
    @test codim(Cone3) == 0
    @test faces(Cone2, 2) isa SubObjectIterator{Cone{T}}
    @test length(faces(Cone2, 2)) == 2
    @test faces(Cone4, 1) isa SubObjectIterator{Cone{T}}
    @test length(faces(Cone4, 1)) == 2
    @test issetequal(faces(Cone2, 2), positive_hull.(Ref(f), [[1 0 0], [0 0 1]], [[0 1 0]]))
    @test _check_im_perm_rows(ray_indices(faces(Cone2, 2)), [[1], [2]])
    @test _check_im_perm_rows(incidence_matrix(faces(Cone2, 2)), [[1], [2]])
    @test _check_im_perm_rows(faces(IncidenceMatrix, Cone2, 2), [[1], [2]])
    @test issetequal(faces(Cone4, 1), positive_hull.(Ref(f), [[0 0 1], [1 0 0]]))
    @test _check_im_perm_rows(ray_indices(faces(Cone4, 1)), [[1], [2]])
    @test _check_im_perm_rows(incidence_matrix(faces(Cone4, 1)), [[1], [2]])
    @test _check_im_perm_rows(faces(IncidenceMatrix, Cone4, 1), [[1], [2]])
    @test _check_im_perm_rows(incidence_matrix(faces(Cone5, 1)), [[1], [2], [3], [4]])
    @test isnothing(faces(Cone2, 1))

    @test f_vector(Cone5) == [4, 4]
    @test f_vector(Cone2) == [0, 2]
    @test lineality_dim(Cone5) == 0
    @test lineality_dim(Cone2) == 1
    @test facet_degrees(Cone5) == fill(2, 4)
    @test facet_degrees(Cone6) == fill(1, 2)
    @test ray_degrees(Cone5) == fill(2, 4)
    @test ray_degrees(Cone6) == fill(1, 2)

    @test n_facets(Cone5) == 4
    @test relative_interior_point(Cone1) == f.([1//2, 1//2])
    @test length(findall(f -> [1, 0, 0] in f, facets(Hyperplane, Cone5))) == 2
    @test length(findall(f -> [1, 0, 0] in f, facets(Halfspace, Cone5))) == 4
  end

  @testset "constructors" begin
    @test cone_from_inequalities(f, [-1 0 0; 0 0 -1]) == Cone2
    @test cone_from_inequalities(f, [-1 0 0; 0 0 -1]; non_redundant=true) == Cone2
    @test cone_from_inequalities(f, facets(Cone4), linear_span(Cone4)) == Cone4
    @test cone_from_inequalities(
      f, facets(Cone4), linear_span(Cone4); non_redundant=true
    ) == Cone4
    @test cone_from_equations(f, [0 1 0]) ==
      cone_from_inequalities(f, Matrix{Int}(undef, 0, 3), linear_span(Cone4))
  end
end

@testset "transform $T" for (f, T) in _prepare_scalar_types()
  pts = [1 0 0 0; 1 1 0 0; 1 1 1 0; 1 0 1 0]
  lin = [0 0 0 1]
  A = [-6 3 1 3;
    -24 7 4 13;
    -31 9 5 17;
    -37 11 6 20]
  # We also put inverse here, since inversion in Julia produces float errors
  # and we don't want to do complicated matrix conversions here.
  invA = [1 2 3 -4;
    1 0 1 -1;
    1 9 0 -6;
    1 1 5 -5]
  C = positive_hull(f, pts, lin)
  Ctarget = positive_hull(f, pts * transpose(A), lin * transpose(A))
  for props in (["RAYS", "LINEALITY_SPACE"],
    ["FACETS", "LINEAR_SPAN"])
    Ccopy = Polymake.polytope.Cone{Oscar._scalar_type_to_polymake(T)}()
    for prop in props
      Polymake.take(Ccopy, prop, Polymake.give(Oscar.pm_object(C), prop))
    end
    Ccopy = Cone{T}(Ccopy, f)
    Ccopyt = transform(Ccopy, A)
    Ccopytt = transform(Ccopyt, invA)
    @test Ccopy == C
    @test Ctarget == Ccopyt
    @test Ccopytt == C
  end
end
