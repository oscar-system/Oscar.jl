@testset "PolyhedralFan{$T}" for (f, T) in _prepare_scalar_types()
  C0 = cube(f, 2)
  @test normal_fan(C0) isa PolyhedralFan{T}
  NFsquare = normal_fan(C0)
  R = [1 0 0; 0 0 1]
  L = [0 1 0]
  Cone4 = positive_hull(f, R)
  Cone5 = positive_hull(f, [1 0 0; 0 1 0])

  @test polyhedral_fan([Cone4, Cone5]) isa PolyhedralFan{T}
  F0 = polyhedral_fan([Cone4, Cone5])
  I3 = [1 0 0; 0 1 0; 0 0 1]
  incidence1 = incidence_matrix([[1, 2], [2, 3]])
  incidence2 = incidence_matrix([[1, 2]])
  @test polyhedral_fan(f, incidence1, I3) isa PolyhedralFan{T}
  F1 = polyhedral_fan(f, incidence1, I3)
  F1NR = polyhedral_fan(f, incidence1, I3; non_redundant=true)
  @test polyhedral_fan(f, incidence1, I3) isa PolyhedralFan{T}
  F2 = polyhedral_fan(f, incidence2, R, L)
  F2NR = polyhedral_fan(f, incidence2, R, L; non_redundant=true)

  @test Polymake.exists(Oscar.pm_object(F2NR), "RAYS")
  @test !Polymake.exists(Oscar.pm_object(F2NR), "INPUT_RAYS")
  @test !Polymake.exists(Oscar.pm_object(F2), "RAYS")
  @test Polymake.exists(Oscar.pm_object(F2), "INPUT_RAYS")

  @testset "core functionality" begin
    if T == QQFieldElem
      @test is_smooth(NFsquare)
    end
    @test rays(NFsquare) isa SubObjectIterator{RayVector{T}}
    @test length(rays(NFsquare)) == 4
    @test issetequal(
      rays(NFsquare), ray_vector.(Ref(f), [[1, 0], [-1, 0], [0, 1], [0, -1]])
    )
    @test vector_matrix(rays(NFsquare)) == _oscar_matrix_from_property(f, rays(NFsquare))
    @test is_regular(NFsquare)
    @test is_complete(NFsquare)
    @test !is_complete(F0)
    @test length(rays(F0)) == 3
    @test n_rays(F1) == 3
    @test dim(F1) == 2
    @test ambient_dim(F1) == 3
    @test n_rays(F2) == 0
    @test lineality_dim(F2) == 1
    RMLF2 = rays_modulo_lineality(F2)
    @test length(RMLF2[:rays_modulo_lineality]) == 2
    @test maximal_cones(F1) isa SubObjectIterator{Cone{T}}
    @test dim.(maximal_cones(F1)) == [2, 2]
    @test _check_im_perm_rows(ray_indices(maximal_cones(F1)), incidence1)
    @test _check_im_perm_rows(incidence_matrix(maximal_cones(F1)), incidence1)
    @test _check_im_perm_rows(maximal_cones(IncidenceMatrix, F1), incidence1)
    @test number_of_maximal_cones(F1) == 2
    @test lineality_space(F2) isa SubObjectIterator{RayVector{T}}
    @test length(lineality_space(F2)) == 1
    @test lineality_space(F2) == [L[:]]
    @test generator_matrix(lineality_space(F2)) == matrix(f, L)
    @test matrix(f, lineality_space(F2)) == matrix(f, L)
    @test cones(F2, 2) isa SubObjectIterator{Cone{T}}
    @test size(cones(F2, 2)) == (2,)
    @test lineality_space(cones(F2, 2)[1]) == [[0, 1, 0]]
    @test rays.(cones(F2, 2)) == [[], []]
    @test _check_im_perm_rows(ray_indices(cones(F1, 2)), incidence1)
    @test _check_im_perm_rows(incidence_matrix(cones(F1, 2)), incidence1)
    @test _check_im_perm_rows(cones(IncidenceMatrix, F1, 2), incidence1)

    # Construct a fan that describes affine 3-space
    A3 = polyhedral_fan(positive_hull(identity_matrix(ZZ, 3)))
    @test length(cones(A3, 1)) == 3
    @test length(cones(A3, 0)) == 1
    @test length(cones(A3, -1)) == 0
    @test cones(A3, 1) isa SubObjectIterator{Cone{QQFieldElem}}
    @test cones(A3, 0) isa SubObjectIterator{Cone{QQFieldElem}}
    @test cones(A3, -1) isa SubObjectIterator{Cone{QQFieldElem}}

    II = ray_indices(maximal_cones(NFsquare))
    NF0 = polyhedral_fan(II, rays(NFsquare))
    @test n_rays(NF0) == 4
    FF0 = face_fan(C0)
    @test n_rays(FF0) == 4
    if T == QQFieldElem
      @test !is_smooth(FF0)
    end
    @test f_vector(NFsquare) == [4, 4]
    @test rays(F1NR) == collect(eachrow(I3))
    @test _check_im_perm_rows(ray_indices(maximal_cones(F1NR)), incidence1)
    @test _check_im_perm_rows(incidence_matrix(maximal_cones(F1NR)), incidence1)
    @test n_rays(F2NR) == 0
    @test lineality_dim(F2NR) == 1
    RMLF2NR = rays_modulo_lineality(F2NR)
    @test length(RMLF2NR[:rays_modulo_lineality]) == 2
    @test RMLF2NR[:rays_modulo_lineality] == collect(eachrow(R))
    @test lineality_space(F2NR) == collect(eachrow(L))
    @test _check_im_perm_rows(ray_indices(maximal_cones(F2NR)), incidence2)
    @test _check_im_perm_rows(incidence_matrix(maximal_cones(F2NR)), incidence2)

    C = positive_hull(f, identity_matrix(ZZ, 0))
    pf = polyhedral_fan(C)
    @test number_of_maximal_cones(pf) == 1
    @test dim(pf) == 0
    pfc = polyhedral_fan(maximal_cones(pf))
    @test number_of_maximal_cones(pfc) == 1

    C = positive_hull(f, vcat(identity_matrix(ZZ, 4), -identity_matrix(ZZ, 4)))
    pf = polyhedral_fan(C)
    @test number_of_maximal_cones(pf) == 1

    @test polyhedral_fan(maximal_cones(F2NR)) isa PolyhedralFan
    # this should just deepcopy the object
    F2NRc = polyhedral_fan(maximal_cones(F2NR))
    @test Polymake.list_properties(Oscar.pm_object(F2NR)) ==
      Polymake.list_properties(Oscar.pm_object(F2NRc))

    # construct from a generic list of cones
    @test polyhedral_fan(cones(F1NR, 1)) isa PolyhedralFan
    @test f_vector(polyhedral_fan(cones(F1NR, 1))) == [3]
    @test f_vector(polyhedral_fan(Cone5)) == [2, 1]

    @test polyhedral_fan(cones(F1NR, 1); non_redundant=true) isa PolyhedralFan
    @test f_vector(polyhedral_fan(cones(F1NR, 1))) == [3]
    @test f_vector(polyhedral_fan(Cone5)) == [2, 1]
  end
end

@testset "minimal supercone" begin
  # Construct a fan that describes affine 3-space
  PF = polyhedral_fan(positive_hull(identity_matrix(ZZ, 3)))
  @test issetequal(minimal_supercone_indices(PF, [1, 1, 1]), [1, 2, 3])
  @test minimal_supercone_coordinates(PF, [1, 1, 1]) == [1, 1, 1]
  @test issetequal(minimal_supercone_indices(PF, [1, 1, 0]), [1, 2])
  @test minimal_supercone_coordinates(PF, [1, 1, 0]) == [1, 1, 0]
  @test issetequal(minimal_supercone_indices(PF, [1, 0, 0]), [1])
  @test minimal_supercone_coordinates(PF, [1, 0, 0]) == [1, 0, 0]
  @test issetequal(minimal_supercone_indices(PF, [0, 0, 0]), Int64[])
  @test minimal_supercone_coordinates(PF, [0, 0, 0]) == [0, 0, 0]
  @test issetequal(minimal_supercone_indices(PF, [2//5, 3, 0]), [1, 2])
  @test minimal_supercone_coordinates(PF, [2//5, 3, 0]) == [2//5, 3, 0]

  # Construct a fan that describes projective 3-space
  PF = normal_fan(Oscar.simplex(3))
  @test issetequal(minimal_supercone_indices(PF, [1, 1, 1]), [1, 2, 3])
  @test minimal_supercone_coordinates(PF, [1, 1, 1]) == [1, 1, 1, 0]
  @test issetequal(minimal_supercone_indices(PF, [1, 1, 0]), [1, 2])
  @test minimal_supercone_coordinates(PF, [1, 1, 0]) == [1, 1, 0, 0]
  @test issetequal(minimal_supercone_indices(PF, [1, 0, 0]), [1])
  @test minimal_supercone_coordinates(PF, [1, 0, 0]) == [1, 0, 0, 0]
  @test issetequal(minimal_supercone_indices(PF, [0, 0, 0]), Int64[])
  @test minimal_supercone_coordinates(PF, [0, 0, 0]) == [0, 0, 0, 0]
  @test issetequal(minimal_supercone_indices(PF, [-1, 1, 1]), [2, 3, 4])
  @test minimal_supercone_coordinates(PF, [-1, 1, 1]) == [0, 2, 2, 1]
  @test issetequal(minimal_supercone_indices(PF, [-1, -1, 1]), [3, 4])
  @test minimal_supercone_coordinates(PF, [-1, -1, 1]) == [0, 0, 2, 1]
  @test issetequal(minimal_supercone_indices(PF, [-2//3, QQFieldElem(-4//3), 1]), [1, 3, 4])
  @test minimal_supercone_coordinates(PF, [-2//3, -4//3, 1]) == [2//3, 0, 7//3, 4//3]

  @test is_minimal_supercone_coordinate_vector(PF, [1, 1, 1, 0])
  @test !is_minimal_supercone_coordinate_vector(PF, [1, 1, 1, 1])
  @test !is_minimal_supercone_coordinate_vector(PF, [1, 1, -1, 0])

  @test standard_coordinates(PF, [1, 2, 3, 0]) == [1, 2, 3]
  @test standard_coordinates(PF, [1, 1, 0, 1]) == [0, 0, -1]
  @test standard_coordinates(PF, [1, 0, 0, 1]) == [0, -1, -1]

  # Construct a nonsimplicial fan
  PF = face_fan(cube(3))
  @test minimal_supercone_coordinates(PF, [1, 0, 0]) == [0, 0, 0, 1//2, 0, 1//2, 0, 0]

  # Construct a fan of the 1/2(1, 1) singularity
  ray_generators = [[2, -1], [0, 1]]
  max_cones = IncidenceMatrix([[1, 2]])
  PF = polyhedral_fan(max_cones, ray_generators)
  @test minimal_supercone_coordinates(PF, [1, 0]) == [1//2, 1//2]
end

@testset "Transform{$T}" for (f, T) in _prepare_scalar_types()
  square = cube(f, 2)
  nf_square = normal_fan(square)
  m_good = matrix(f, [1 0; 0 1])
  m_bad = matrix(f, [1 1])
  nf_transformed = transform(nf_square, m_good)
  @test nf_transformed isa PolyhedralFan{T}
  @test n_rays(nf_transformed) == 4
  @test number_of_maximal_cones(nf_transformed) == 4
  @test_throws ErrorException transform(nf_square, m_bad)
  nf_unverified = transform(nf_square, m_good; check=false)
  @test nf_unverified isa PolyhedralFan{T}
  @test n_rays(nf_unverified) == 4
  @test number_of_maximal_cones(nf_unverified) == 4

  cdeg = convex_hull(f, [1 0 0; 0 0 0])
  nflin = normal_fan(cdeg)
  mm = matrix(f, [1 0 0; 0 1 0; 0 0 1])
  nflin_transformed = transform(nflin, mm)
  @test nflin_transformed isa PolyhedralFan{T}
  @test n_rays(nflin_transformed) == 0
  @test number_of_maximal_cones(nflin_transformed) == 2
  @test lineality_dim(nflin_transformed) == 2
end

@testset "Star Subdivision" begin
  f = polyhedral_fan(incidence_matrix([[1, 2, 3, 4]]), [1 0 0; 1 1 0; 1 1 1; 1 0 1])
  @test is_pure(f)
  @test is_fulldimensional(f)
  v0 = [1; 0; 0]
  v1 = [2; 1; 0]
  v2 = [2; 1; 1]
  sf0 = star_subdivision(f, v0)
  sf1 = star_subdivision(f, v1)
  sf2 = star_subdivision(f, v2)
  @test number_of_maximal_cones(sf0) == 2
  @test number_of_maximal_cones(sf1) == 3
  @test number_of_maximal_cones(sf2) == 4

  ff = polyhedral_fan(incidence_matrix([[1], [2, 3]]), [1 0 0; -1 0 0; 0 1 0])
  @test !is_pure(ff)
  @test !is_fulldimensional(ff)
  w0 = [1; 0; 0]
  w1 = [-1; 1; 0]
  sff0 = star_subdivision(ff, w0)
  sff1 = star_subdivision(ff, w1)
  @test number_of_maximal_cones(sff0) == 2
  @test number_of_maximal_cones(sff1) == 3

  f = normal_fan(cube(2))
  fMinus = -f
  @test n_rays(fMinus) == n_rays(f)
  @test n_maximal_cones(fMinus) == n_maximal_cones(f)

  fMinus = -1*f
  @test n_rays(fMinus) == n_rays(f)
  @test n_maximal_cones(fMinus) == n_maximal_cones(f)
end

@testset "Defining polynomial of hyperplane arrangement from matrix" begin
  A = identity_matrix(QQ, 3)
  L = arrangement_polynomial(A)
  x = gens(parent(L))
  @test L == x[1] * x[2] * x[3]
  R, y = polynomial_ring(QQ, [:x, :y, :z])
  LL = arrangement_polynomial(R, A)
  @test LL == y[1] * y[2] * y[3]
  A = identity_matrix(QQ, 2)
  AA = identity_matrix(GF(3), 3)
  @test_throws ArgumentError arrangement_polynomial(R, A)
  @test_throws ArgumentError arrangement_polynomial(R, AA)
  A = matrix(QQ, [1 1 0])
  Avv = [A[i, :] for i in 1:1]
  AM = Matrix(A)
  apAvv = arrangement_polynomial(Avv)
  RAvv = parent(apAvv)
  apAM = arrangement_polynomial(RAvv, AM)
  @test nvars(RAvv) == 3
  @test apAvv == apAM
  @test nvars(parent(arrangement_polynomial(AM))) == 3
  @test nvars(parent(arrangement_polynomial(RAvv, Avv))) == 3
  @test arrangement_polynomial(RAvv, Avv) == apAvv
end
