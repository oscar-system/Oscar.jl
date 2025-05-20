@testset "Scalar types" begin

  # small helper for quantitative testing
  # returns the scalar product of the difference of the vertices of an edge x (which is encoded as a Polyhedron),
  # i.e. the square of the edge length
  function _edge_length_for_test(x::Polyhedron)
    v = vertices(x)
    s = v[1] - v[2]
    return sum([s[i]^2 for i in 1:3])
  end

  # this example relies on the same source as polymake's johnson_solid(84):
  # https://de.wikipedia.org/wiki/Trigondodekaeder
  # as of now, johnson solid J84 is realized as a Polyhedron over floats within polymake
  Qx, x = QQ[:x]
  K, r = number_field(x^3 - 3x^2 - 4x + 8, "r")
  Ky, y = K[:y]
  L, lq = number_field(y^2 - (2 - r^2)//2, "q")
  Lz, z = L[:z]
  emb = only(filter(x -> is_positive(lq, x), real_embeddings(L)))
  E, q = Hecke.embedded_field(L, emb)
  pq1 = only(filter(x -> x > 0, E.(roots(z^2 - (3 + 2r - r^2)))))
  pq2 = only(filter(x -> x > 0, E.(roots(z^2 - (3 - r^2)))))
  p = (pq1 + pq2)//2
  r = E(r)

  e = one(QQBarFieldElem)
  s, c = sinpi(2 * e / 5), cospi(2 * e / 5)
  pts = [[e, e], [s, c], [c, s]]

  v = [0 1 p; 0 -1 p; r 0 q; -r 0 q; 0 r -q; 0 -r -q; 1 0 -p; -1 0 -p]
  V = matrix(E, v)
  sd = convex_hull(E, v; non_redundant=true)

  @test sd isa Polyhedron{
    Hecke.EmbeddedNumFieldElem{Hecke.RelSimpleNumFieldElem{AbsSimpleNumFieldElem}}
  }
  @test point_matrix(vertices(sd)) == V
  # volume formula from source
  @test volume(sd) == 8r//12 * (6pq2 + 2 * r * q)
  # the snub disphenoid consists of 12 triangles
  @test n_vertices.(faces(sd, 2)) == repeat([3], 12)
  # and here the 18 edges are of length 2
  @test _edge_length_for_test.(faces(sd, 1)) == repeat([4], 18)
  # scaling the Polyhedron by 3 yields edge lengths of 6
  @test _edge_length_for_test.(faces(3 * sd, 1)) == repeat([36], 18)
  # there are 11 lattice points
  @test length(lattice_points(sd)) == 11

  let pc = polyhedral_complex(
      E, incidence_matrix(facets(sd)), vertices(sd); non_redundant=true
    )
    @test issetequal(maximal_polyhedra(pc), faces(sd, 2))
  end
  let c = convex_hull(E, permutedims([0]), permutedims([r]))
    ms = product(sd, c)
    @test recession_cone(ms) == positive_hull(E, [0 0 0 1])
  end
  nf = normal_fan(sd)
  nfc = polyhedral_fan(E, maximal_cones(IncidenceMatrix, nf), rays(nf))
  @test is_regular(nfc)

  @test convex_hull(pts) isa Polyhedron{QQBarFieldElem}
  qp = convex_hull(pts)

  @test number_of_vertices(qp) == 3
  @test number_of_facets(qp) == 3

  @test length(lattice_points(qp)) == 1

  @testset "Scalar detection" begin
    let j = johnson_solid(12)
      @test j isa Polyhedron{QQFieldElem}
      s = vertices(j)[1][1]
      @test s isa QQFieldElem
    end
    let j = polyhedron(Polymake.polytope.johnson_solid(64))
      @test j isa Polyhedron{Float64}
      s = vertices(j)[1][1]
      @test s isa Float64
    end
    let a = archimedean_solid("truncated_cube")
      @test a isa Polyhedron{Hecke.EmbeddedNumFieldElem{AbsSimpleNumFieldElem}}
      s = vertices(a)[1][1]
      @test s isa Hecke.EmbeddedNumFieldElem{AbsSimpleNumFieldElem}
      f = coefficient_field(a)
      F = number_field(f)
      isq = Hecke.is_quadratic_type(F)
      @test isq[1]
      @test isq[2] == 2
    end
    let j = johnson_solid(1)
      jj = polyhedron(
        Polymake.polytope.Polytope{Polymake.OscarNumber}(;
          POINTS=Oscar.pm_object(j).VERTICES
        ),
      )
      @test number_field(coefficient_field(j)) == number_field(coefficient_field(jj))
    end
    let ng = n_gon(5)
      (A, b) = halfspace_matrix_pair(facets(ng))
      @test typeof(polyhedron(A, b)) == typeof(ng)
      @test coefficient_field(polyhedron(A, b)) == coefficient_field((ng))
    end
  end

  @testset "Cross scalar operations" begin
    NF, sr2 = quadratic_field(2)
    ENF, sre2 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
    T = elem_type(ENF)

    c = cube(ENF, 3, -sre2, sre2)
    d = cube(3, -5//3, 1//2)

    mat = [1 2 3; 1 4 5]
    e = positive_hull(ENF, mat; non_redundant=true)
    f = positive_hull(mat; non_redundant=true)

    # Polyhedron
    for z in (4, QQ(4), ENF(4))
      @test pyramid(c, z) isa Polyhedron{T}
      p = pyramid(c, z)
      @test f_vector(p) == [9, 20, 18, 7]
      @test volume(p) == volume(c)
    end
    for (z, U) in ((4, QQFieldElem), (QQ(4), QQFieldElem), (ENF(4), T))
      @test pyramid(d, z) isa Polyhedron{U}
      p = pyramid(d, z)
      @test f_vector(p) == [9, 20, 18, 7]
      @test volume(p) == volume(d)
    end

    for z in (2, QQ(2), ENF(2))
      @test bipyramid(c, z) isa Polyhedron{T}
      p = bipyramid(c, z)
      @test f_vector(p) == [10, 28, 30, 12]
      @test volume(p) == volume(c)
      for z_prime in (-2, QQ(-2), ENF(-2))
        @test bipyramid(c, z, z_prime) isa Polyhedron{T}
        p = bipyramid(c, z, z_prime)
        @test f_vector(p) == [10, 28, 30, 12]
        @test volume(p) == volume(c)
      end
    end
    for (z, U) in ((2, QQFieldElem), (QQ(2), QQFieldElem), (ENF(2), T))
      @test bipyramid(d, z) isa Polyhedron{U}
      p = bipyramid(d, z)
      @test f_vector(p) == [10, 28, 30, 12]
      @test volume(p) == volume(d)
      for (z_prime, V) in ((-2, QQFieldElem), (QQ(-2), QQFieldElem), (ENF(-2), T))
        @test bipyramid(d, z, z_prime) isa Polyhedron{U == T || V == T ? T : U}
        p = bipyramid(d, z, z_prime)
        @test f_vector(p) == [10, 28, 30, 12]
        @test volume(p) == volume(d)
      end
    end

    @test intersect(c, d) isa Polyhedron{T}
    let p = intersect(c, d)
      @test f_vector(p) == f_vector(c)
      @test volume(p) == 11 * sre2//4 + 25//8
    end

    @test minkowski_sum(c, d) isa Polyhedron{T}
    let p = minkowski_sum(c, d)
      @test f_vector(p) == f_vector(c)
      @test volume(p) == 265 * sre2//6 + 13429//216
      @test c + d == p
    end

    for (k, U) in ((3, QQFieldElem), (QQ(3), QQFieldElem), (ENF(3), T))
      @test k * c isa Polyhedron{T}
      let p = k * c
        @test f_vector(p) == f_vector(c)
        @test volume(p) == 27 * volume(c)
        @test c * k == p
      end
      @test k * d isa Polyhedron{U}
      let p = k * d
        @test f_vector(p) == f_vector(d)
        @test volume(p) == 27 * volume(d)
        @test d * k == p
      end
    end

    @test convex_hull(c, d) isa Polyhedron{T}
    let p = convex_hull(c, d)
      @test f_vector(p) == [14, 24, 12]
      @test volume(p) == 379 * sre2//36 + 1349//108
    end

    cc = positive_hull(ENF, Oscar.homogenized_matrix(vertices(c), 1))
    dc = positive_hull(Oscar.homogenized_matrix(vertices(d), 1))
    @test intersect(cc, dc) isa Cone{T}
    let p = intersect(cc, dc)
      @test f_vector(p) == f_vector(c)
      @test coefficient_field(p) == ENF
    end

    tm = [0 0 1; 0 1 0; 1 0 0]
    # Cone
    for (m, F) in (
      (tm, QQ),
      (Rational{BigInt}.(tm), QQ),
      (QQ.(tm), QQ),
      (matrix(QQ, tm), QQ),
      (ENF.(tm), ENF),
      (matrix(ENF, tm), ENF),
    )
      U = elem_type(F)
      @test transform(e, m) isa Cone{T}
      let et = transform(e, m)
        @test f_vector(et) == f_vector(e)
        @test vector_matrix(rays(et)) == matrix(ENF, [1 2//3 1//3; 1 4//5 1//5])
      end
      @test transform(f, m) isa Cone{U}
      let ft = transform(f, m)
        @test f_vector(ft) == f_vector(f)
        @test vector_matrix(rays(ft)) == matrix(F, [1 2//3 1//3; 1 4//5 1//5])
      end
    end
  end

  @testset "QuadraticExtension-templated sub-objects" begin
    j = johnson_solid(2)

    let f = facets(Polyhedron, j)
      fj = normal_vector.(facets(j))
      fj = [fj; -fj]
      for i in 1:n_facets(j)
        @test normal_vector(affine_hull(f[i])[]) in fj
      end
      @test halfspace_matrix_pair(f) isa NamedTuple{
        (:A, :b),
        Tuple{
          AbstractAlgebra.Generic.MatSpaceElem{EmbeddedNumFieldElem{AbsSimpleNumFieldElem}},
          Vector{EmbeddedNumFieldElem{AbsSimpleNumFieldElem}},
        },
      }
      g = halfspace_matrix_pair(f)
      @test affine_halfspace(coefficient_field(j), g.A[1, :], g.b[1]) in facets(j)
    end
    for n in (1, 2) # faces which are facets have a different access function
      let f = faces(j, n)
        for i in 1:Int(f_vector(j)[n + 1])
          @test issubset(vertices(f[i]), vertices(j))
        end
      end
    end

    k = face_fan(j)
    let f = cones(k, 3)
      for i in 1:Int(f_vector(k)[3])
        @test issubset(rays(f[i]), rays(k))
      end
      l = f[1]
    end
  end
end
