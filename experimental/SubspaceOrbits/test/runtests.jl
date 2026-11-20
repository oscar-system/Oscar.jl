@testset "SubspaceOrbits" begin
  SO = Oscar.SubspaceOrbits

  # Orbits on the labels 1:D computed with the reference implementation
  # and a plain orbit algorithm: `orbit_of[l]` is the number of the orbit
  # of `l`, and `lengths` are the orbit lengths.
  function naive_label_orbits(C::SO.LabelCodec, gens::Vector{Matrix{Int}})
    orbit_of = zeros(Int, C.D)
    lengths = Int[]
    for l in 1:C.D
      orbit_of[l] == 0 || continue
      push!(lengths, 1)
      orbit_of[l] = length(lengths)
      queue = [l]
      while !isempty(queue)
        u = pop!(queue)
        for g in gens
          v = SO.reference_image(C, g, u)
          orbit_of[v] == 0 || continue
          orbit_of[v] = length(lengths)
          lengths[end] += 1
          push!(queue, v)
        end
      end
    end
    return orbit_of, lengths
  end

  # The label of a subspace of the natural module.
  function label_of_subspace(C::SO.LabelCodec, U)
    M = SO.int_matrix(rref(matrix(U.map))[2])
    return SO.matrix_to_label(C, M)
  end

  # Compare both user functions with the naive orbit computation.
  function check_orbits(G, k)
    d = degree(G)
    p = Int(order(base_ring(G)))
    C = SO.LabelCodec(p, d, k)
    orbit_of, lengths = naive_label_orbits(C, [SO.int_matrix(matrix(g)) for g in gens(G)])

    res = subspace_orbit_representatives_and_lengths(G, k)
    ids = [orbit_of[label_of_subspace(C, U)] for (U, n) in res]
    @test allunique(ids)
    @test length(ids) == length(lengths)
    @test all(lengths[ids[i]] == res[i][2] for i in 1:length(res))

    res2 = subspace_orbit_representatives_and_stabilizers(G, k)
    @test [U for (U, n) in res] == [U for (U, S) in res2]
    for i in 1:length(res2)
      U, S = res2[i]
      l = label_of_subspace(C, U)
      @test order(S) * res[i][2] == order(G)
      # the order of `S` may have been set, not computed
      @test order(sub(G, gens(S))[1]) == order(S)
      @test all(SO.reference_image(C, SO.int_matrix(matrix(s)), l) == l for s in gens(S))
      @test is_subgroup(S, G)[1]
    end
  end

  # Run the orbit algorithm on a small case with what large cases select:
  # the radical also for lengths, dropped blocks, sampled edges, 64 bit labels.
  function check_large_case_branches(G, k)
    C = SO.LabelCodec(Int(order(base_ring(G))), degree(G), k)
    orbit_of, lengths = naive_label_orbits(C, [SO.int_matrix(matrix(g)) for g in gens(G)])
    expected = [(findfirst(==(i), orbit_of), lengths[i]) for i in 1:length(lengths)]
    for stabilizers in (false, true), L in (UInt32, UInt64)
      A = SO.SubspaceAction(G, k; stabilizers, min_images_saved=0)
      got = Tuple{Int,Int}[]
      SO.label_orbits(A.action, A.relorders; track=stabilizers, trim_threshold=4,
                      label_type=L) do state
        push!(got, (state.rep, state.length))
        stabilizers || return nothing
        S = SO.point_stabilizer(A, state; max_shuffled_edges=0)
        @test order(sub(G, gens(S))[1]) * state.length == order(G)
      end
      @test got == expected
    end
  end

  @testset "labels for p = $p, d = $d, k = $k" for (p, d, k) in
      [(2, 4, 2), (2, 5, 3), (2, 6, 2), (3, 3, 2), (3, 4, 2), (5, 3, 1), (5, 3, 2), (7, 2, 1)]
    C = SO.LabelCodec(p, d, k)
    @test C.D == Int(q_binomial(d, k, ZZ(p)))
    mats = [SO.label_to_matrix(C, l) for l in 1:C.D]
    @test allunique(mats)
    @test all(M -> SO.rref_mod_p!(copy(M), p) == M, mats)
    @test all(l -> SO.matrix_to_label(C, mats[l]) == l, 1:C.D)
  end

  # The last three cases exercise F_2 rows in several chunks, sums of rows
  # that need intermediate reductions, and 128 bit rows.
  @testset "label action for p = $p, d = $d, k = $k" for (p, d, k) in
      [(2, 6, 3), (2, 10, 4), (3, 5, 2), (3, 8, 3), (5, 4, 2), (7, 4, 2), (251, 3, 1),
       (2, 14, 3), (3, 14, 1), (251, 5, 1)]
    C = SO.LabelCodec(p, d, k)
    F = GF(p)
    gens = [SO.int_matrix(matrix(rand(GL(d, F)))) for i in 1:3]
    K = SO.LabelAction(C, gens)
    labels = C.D <= 1000 ? (1:C.D) : rand(1:C.D, 1000)
    @test all(SO.image(K, i, l) == SO.reference_image(C, gens[i], l)
              for i in 1:3, l in labels)
  end

  @testset "orbits" begin
    # solvable groups
    G = sylow_subgroup(GL(4, 3), 3)[1]
    for k in 1:3
      check_orbits(G, k)
    end
    G = sylow_subgroup(GL(5, 2), 2)[1]
    for k in 1:4
      check_orbits(G, k)
    end
    G = matrix_group(permutation_matrix.(Ref(GF(5)), gens(symmetric_group(4))))
    for k in 1:3
      check_orbits(G, k)
    end

    # non-solvable groups, with trivial and with nontrivial solvable radical
    G = matrix_group(permutation_matrix.(Ref(GF(2)), gens(symmetric_group(5))))
    for k in 1:4
      check_orbits(G, k)
    end
    check_orbits(GL(4, 2), 2)
    check_orbits(GL(3, 5), 2)
    check_orbits(GL(3, 3), 1)
    G = stabilizer(GL(5, 2), [GF(2)(x) for x in [1, 0, 0, 0, 0]], *)[1]
    for k in 1:4
      check_orbits(G, k)
    end

    # redundant generators, trivial group
    c = permutation_matrix(GF(2), cperm(1:6))
    check_orbits(matrix_group([c, c^2, c^3, c^4, c^5]), 2)
    G = matrix_group([identity_matrix(GF(2), 4)])
    @test length(subspace_orbit_representatives_and_lengths(G, 2)) == 35
    check_orbits(G, 2)
  end

  @testset "branches for large cases" begin
    check_large_case_branches(GL(4, 2), 2)
    check_large_case_branches(GL(3, 5), 1)
    check_large_case_branches(sylow_subgroup(GL(5, 2), 2)[1], 2)
    G = stabilizer(GL(5, 2), [GF(2)(x) for x in [1, 0, 0, 0, 0]], *)[1]
    check_large_case_branches(G, 2)
  end

  @testset "edge cases" begin
    G = GL(3, 3)
    for k in [0, 3]
      res = subspace_orbit_representatives_and_lengths(G, k)
      @test length(res) == 1 && dim(res[1][1]) == k && res[1][2] == 1
      res = subspace_orbit_representatives_and_stabilizers(G, k)
      @test length(res) == 1 && dim(res[1][1]) == k && res[1][2] == G
    end
    @test_throws ArgumentError subspace_orbit_representatives_and_lengths(G, 4)
    for k in 0:2
      @test_throws ArgumentError subspace_orbit_representatives_and_lengths(GL(2, 4), k)
      @test_throws ArgumentError subspace_orbit_representatives_and_stabilizers(GL(2, 4), k)
    end
    @test_throws ArgumentError subspace_orbit_representatives_and_lengths(GL(30, 5), 15)
    @test_throws ArgumentError subspace_orbit_representatives_and_lengths(GL(2, GF(65537)), 1)

    # all representations of a prime field
    for F in [GF(3), Nemo.Native.GF(3), Nemo.Native.GF(ZZ(3)),
              Nemo.Native.finite_field(3, 1)[1], Nemo.Native.finite_field(ZZ(3), 1)[1]]
      G = GL(2, F)
      @test [n for (U, n) in subspace_orbit_representatives_and_lengths(G, 1)] == [4]
      res = subspace_orbit_representatives_and_stabilizers(G, 1)
      @test [index(G, S) for (U, S) in res] == [4]
    end

    # no group theory for large modules
    res = subspace_orbit_representatives_and_lengths(GL(2, GF(65521)), 1)
    @test [n for (U, n) in res] == [65522]
    @test_throws ArgumentError subspace_orbit_representatives_and_stabilizers(GL(2, GF(4099)), 1)

    # a prime field that GAP cannot be asked about
    F = Oscar.AbstractAlgebra.GF(3)
    G = matrix_group([matrix(F, [0 1; 1 0])])
    @test sort([n for (U, n) in subspace_orbit_representatives_and_lengths(G, 1)]) == [1, 1, 2]
    @test_throws ArgumentError subspace_orbit_representatives_and_stabilizers(G, 1)

    # too many subspaces are rejected without counting them
    err = @test_throws ArgumentError SO.LabelCodec(2, 2000, 1000)
    @test length(err.value.msg) < 100

    # labels must fit into the type they are stored in
    K = SO.LabelAction(SO.LabelCodec(2, 9, 1), Matrix{Int}[])
    @test_throws ArgumentError SO.label_orbits(identity, K, Int[]; track=true, label_type=UInt8)

    # no group theory where it cannot pay off
    G = matrix_group(permutation_matrix.(Ref(GF(2)), gens(symmetric_group(5))))
    @test length(subspace_orbit_representatives_and_lengths(G, 2)) == 10
    @test !has_order(G)
  end

  @testset "consistency with orbit_representatives_and_stabilizers" begin
    G = orthogonal_group(1, 4, GF(3))
    for k in 1:3
      res = subspace_orbit_representatives_and_stabilizers(G, k)
      ref = orbit_representatives_and_stabilizers(G, k)
      @test typeof(res) == typeof(ref)
      @test sort([order(S) for (U, S) in res]) == sort([order(S) for (U, S) in ref])
    end

    # a group that is not conjugate to its transpose: it fixes a vector,
    # the transposed group a hyperplane
    G = stabilizer(GL(4, 2), [GF(2)(x) for x in [1, 0, 0, 0]], *)[1]
    @test sort([n for (U, n) in subspace_orbit_representatives_and_lengths(G, 1)]) == [1, 14]
    @test sort([n for (U, n) in subspace_orbit_representatives_and_lengths(G, 3)]) == [7, 8]
    for k in 1:3
      res = subspace_orbit_representatives_and_stabilizers(G, k)
      ref = orbit_representatives_and_stabilizers(G, k)
      @test sort([order(S) for (U, S) in res]) == sort([order(S) for (U, S) in ref])
      for (U, S) in res
        m = rref(matrix(U.map))[2]
        @test S == stabilizer(G, m, on_echelon_form_mats)[1]
        @test length(orbit(G, on_echelon_form_mats, m)) == index(G, S)
      end
    end
  end

  # The label of a subspace is the position of its basis in
  # `bases_of_subspaces`, and the action on labels is `on_echelon_form_mats`.
  @testset "consistency with bases_of_subspaces for p = $p, d = $d, k = $k" for (p, d, k) in
      [(2, 4, 2), (3, 4, 1), (3, 4, 3), (5, 3, 2)]
    F = GF(p)
    G = GL(d, F)
    C = SO.LabelCodec(p, d, k)
    mats = collect(Oscar.bases_of_subspaces(vector_space(F, d), k))
    @test [SO.int_matrix(m) for m in mats] == [SO.label_to_matrix(C, l) for l in 1:C.D]
    K = SO.LabelAction(C, [SO.int_matrix(matrix(g)) for g in gens(G)])
    @test all(mats[SO.image(K, i, l)] == on_echelon_form_mats(mats[l], gen(G, i))
              for i in 1:ngens(G), l in 1:C.D)
  end
end
