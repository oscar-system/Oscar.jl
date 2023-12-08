###############################################################################
#
#  General interface for enumeration of genera for definite ZZLat
#
###############################################################################

#=====TODO
- Implement a serialisation process to save and upload partial results;
- Include the possibility to use isometry enumeration;
======#

###############################################################################
#
#  Neighbours
#
###############################################################################

function neighbour(L::ZZLat, v::ZZMatrix, p::ZZRingElem)
  M = gram_matrix(L)*v
  r, K = left_kernel(matrix(GF(p), rank(L), 1, collect(M)))
  LL = lattice_in_same_ambient_space(L, lift(view(K, 1:r, :))*basis_matrix(L)) + lattice_in_same_ambient_space(L, 1//p*(transpose(v)*basis_matrix(L))) + p*L
  return LL
end

function make_admissible!(w::ZZMatrix, form::ZZMatrix, p::ZZRingElem, K::FinField)
  m = (transpose(w)*form*w)[1]
  a = mod(m, p^2)
  iszero(a) && return nothing
  a = K(div(a, p))
  v = 2*form*w
  v = map_entries(b -> K(b)//a, v)
  v = reduce(vcat, dense_matrix_type(K)[identity_matrix(K, 1), v])
  _, L = left_kernel(v)
  @hassert :ZZLatWithIsom 1 !iszero(view(L, :, 1))
  j = findfirst(j -> !iszero(L[j, 1]), 1:nrows(L))
  @hassert :ZZLatWithIsom 1 !isnothing(j)
  L = map_entries(b -> b//L[j, 1], L)
  v = lift(L[j, 2:ncols(L)])
  for i in 1:nrows(w)
    w[i, 1] += p*v[i]
  end
  @hassert :ZZLatWithIsom 1 iszero(mod(transpose(w)*form*w, p^2))
  return nothing
end

function _neighbours_definite_ex(L::ZZLat, p::ZZRingElem; callback::Function,
                                                          inv_dict::Dict,
                                                          _invariants::Function,
                                                          save_partial::Bool = false,
                                                          save_path::Union{Nothing, IO, String} = nothing,
                                                          use_mass::Bool = true,
                                                          missing_mass::Union{Nothing, Base.RefValue{QQFieldElem}} = nothing)
  K = GF(p)
  Tp = torsion_quadratic_module(L, p*L)
  _, jp = radical_bilinear(Tp)
  form = map_entries(ZZ, gram_matrix(L))

  if use_mass
    __mass = missing_mass[]
  end

  P = Hecke.enumerate_lines(K, rank(L))

  result = typeof(L)[]

  @vprintln :ZZLatWithIsom 1 "$(length(P)) lines to try"

  for x in P

    w = lift(matrix(x))
    a = Tp(abelian_group(Tp)(w))
    iszero(quadratic_product(a)) || continue
    has_preimage(jp, a)[1] && continue

    make_admissible!(w, form, p, K)
    LL = lll(neighbour(L, w, p))

    @hassert :ZZLatWithIsom 1 is_locally_isometric(LL, L, p)

    keep = callback(LL)
    !keep && continue
    @vprintln :ZZLatWithIsom 1 "Keep an isometry class"
    invLL = _invariants(LL)
    if haskey(inv_dict, invLL)
      push!(inv_dict[invLL], LL)
    else
      inv_dict[invLL] = ZZLat[LL]
    end
    push!(result, LL)

    if use_mass
      s = automorphism_group_order(LL)
      sub!(__mass, __mass, 1//s)
      println(__mass)
      is_zero(__mass) && return result
    end
  end
  return result
end

function _neighbours_definite_orbit(L::ZZLat, p::ZZRingElem; callback::Function,
                                                             inv_dict::Dict,
                                                             _invariants::Function,
                                                             save_partial::Bool = false,
                                                             save_path::Union{Nothing, IO, String} = nothing,
                                                             use_mass::Bool = true,
                                                             missing_mass::Union{Nothing, Base.RefValue{QQFieldElem}} = nothing)
  K = GF(p)
  q = GAP.Obj(p)
  form = map_entries(ZZ, gram_matrix(L))
  B = basis_matrix(L)
  Tp = torsion_quadratic_module(L, p*L)
  _, jp = radical_bilinear(Tp)

  if use_mass
    __mass = missing_mass[]
  end

  G = isometry_group(L)
  gensp = dense_matrix_type(K)[map_entries(K, solve_left(B, B*matrix(g))) for g in gens(G)]
  filter!(!is_diagonal, gensp)
  Gp = matrix_group(gensp)
  orbs = GAP.Globals.Orbits(Gp.X, GAP.Globals.Subspaces(GAP.Globals.GF(q)^(rank(L)), 1))
  LO = [[K(x) for x in GAP.Globals.BasisVectors(GAP.Globals.Basis(orb[1]))[1]] for orb in orbs]::Vector{Vector{elem_type(K)}}

  result = typeof(L)[]

  @vprintln :ZZLatWithIsom 1 "$(length(LO)) orbits of lines to try"

  for x in LO

    w = lift(matrix(x))
    a = Tp(abelian_group(Tp)(w))
    iszero(quadratic_product(a)) || continue
    has_preimage(jp, a)[1] && continue

    make_admissible!(w, form, p, K)
    LL = lll(neighbour(L, w, p))

    @hassert :ZZLatWithIsom 1 is_locally_isometric(LL, L, p)

    keep = callback(LL)
    !keep && continue
    @vprintln :ZZLatWithIsom 1 "Keep an isometry class"
    invLL = _invariants(LL)
    if haskey(inv_dict, invLL)
      push!(inv_dict[invLL], LL)
    else
      inv_dict[invLL] = ZZLat[LL]
    end
    push!(result, LL)

    if use_mass
      s = automorphism_group_order(LL)
      sub!(__mass, __mass, 1//s)
      println(__mass)
      is_zero(__mass) && return result
    end
  end
  return result
end

function _neighbours_definite_rand(L::ZZLat, p::ZZRingElem; rand_neigh::Union{Nothing, Int} = nothing,
                                                            callback::Function,
                                                            inv_dict::Dict,
                                                            _invariants::Function,
                                                            save_partial::Bool = false,
                                                            save_path::Union{Nothing, IO, String} = nothing,
                                                            use_mass::Bool = true,
                                                            missing_mass::Union{Nothing, Base.RefValue{QQFieldElem}} = nothing)
  K = GF(p)
  form = map_entries(ZZ, gram_matrix(L))
  Tp = torsion_quadratic_module(L, p*L)
  _, jp = radical_bilinear(Tp)

  if use_mass
    __mass = missing_mass[]
  end

  P = Hecke.enumerate_lines(K, rank(L))

  result = typeof(L)[]

  maxlines = isnothing(rand_neigh) ? min(100, length(P)) : min(rand_neigh, length(P))

  @vprintln :ZZLatWithIsom 1 "Try $(maxlines) random lines"

  for i in 1:maxlines

    w = lift(matrix(rand(P)))
    a = Tp(abelian_group(Tp)(w))
    iszero(quadratic_product(a)) || continue
    has_preimage(jp, a)[1] && continue

    make_admissible!(w, form, p, K)
    LL = lll(neighbour(L, w, p))
    println("bla")
    @hassert :ZZLatWithIsom 1 is_locally_isometric(LL, L, p)

    keep = callback(LL)
    !keep && continue
    @vprintln :ZZLatWithIsom 1 "Keep an isometry class"
    invLL = _invariants(LL)
    if haskey(inv_dict, invLL)
      push!(inv_dict[invLL], LL)
    else
      inv_dict[invLL] = ZZLat[LL]
    end
    push!(result, LL)

    if use_mass
      s = automorphism_group_order(LL)
      sub!(__mass, __mass, 1//s)
      println(__mass)
      is_zero(__mass) && return result
    end
  end
  return result
end

#======User interface
Input:
- known -> finite list of known isometry classes (always non-empty by starting from a single lattice)
- alg_type -> how to enumerate neighbours: all of them (:exhaustive), orbits of them (:orbit), a random number (:random)
Optional:
- rand_neigh -> for random enumeration, how many randome neighbours are computed
- distinct -> if the lattices in "known" are known to be pairwise non-isometric
- invariant_func -> functions to compute isometry invariants for comparing lattices
- save_partial -> whether one wants to save iteratively new isometry classes (for instance for large genera)
- save_path -> in the case "save_partial" is true, where to save the lattices
- use_mass -> whether to use the mass formula as termination condition
- missing_mass -> if "use_mass" and "distinct" are true, and the partial mass of "known" is known, mention what is the part of the mass missing
======#

function _unique_iso_class!(L::Vector{ZZLat})
  isempty(A) && return A
  idxs = eachindex(A)
  y = first(A)
  T = NTuple{2, Any}
  it = iterate(idxs, (iterate(idxs)::T)[2])
  count = 1
  for x in Iterators.drop(A, 1)
    if !is_isometric(x, y)
      it = it::T
      y = A[it[1]] = x
      count += 1
      it = iterate(idxs, it[2])
    end
  end
  resize!(A, count)::typeof(A)
end

function default_func(L::ZZLat)
  m = minimum(L)
  rlr, _ = root_lattice_recognition(L)
  sort!(rlr; lt=(a,b) -> a[1] < b[1] || a[1] == b[1] && a[2] <= b[2])
  kn = kissing_number(L)::Int
  igo = automorphism_group_order(L)::ZZRingElem
  return (m, rlr, kn, igo)
end

function enumerate_definite_genus(
    known::Vector{ZZLat},
    alg_type::Symbol = :exhaustive;
    rand_neigh::Union{Int, Nothing} = nothing,
    distinct::Bool = true,
    invariant_func::Function = default_func,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true,
    _missing_mass::Union{QQFieldElem, Nothing} = nothing
  )
  @req !is_empty(known) "Should know at least one lattice in the genus"
  @req all(LL -> genus(LL) == genus(known[1]), known) "Known lattices must be in the same genus"

  res = copy(known)
  !distinct && _unique_iso_class!(res)

  L, itk = Iterators.peel(res)
  inv_lat = invariant_func(L)
  inv_dict = Dict{typeof(inv_lat), Vector{ZZLat}}(inv_lat => ZZLat[L])
  for N in itk
    inv_lat = invariant_func(N)
    if haskey(inv_dict, inv_lat)
      push!(inv_dict[inv_lat], N)
    else
      inv_dict[inv_lat] = ZZLat[N]
    end
  end

  function _invariants(M::ZZLat)
    for (I, Z) in keys(inv_dict)
      M in Z && return I
    end
    return invariant_func(M)
  end

  callback = function(M::ZZLat)
    any(isequal(M), known) && return false
    invM = _invariants(M)
    !haskey(inv_dict, invM) && return true
    keep = all(N -> !is_isometric(N, M), inv_dict[invM])
    return keep
  end

  if use_mass
    _mass = mass(L)
    if isnothing(_missing_mass)
      found = sum(1//automorphism_group_order(M) for M in res)
      missing_mass = Ref{QQFieldElem}(_mass-found)
    else
      missing_mass = Ref{QQFieldElem}(_missing_mass)
      @hassert :ZZLatWithIsom 1 missing_mass[] <= _mass
    end
  end

  ps = primes(genus(L))
  p = ZZ(3)
  if p in ps
    p = next_prime(p)
  end

  tbv = trues(length(res))
  println(missing_mass[])
  while any(tbv)
    i = findfirst(tbv)
    tbv[i] = false
    if alg_type == :exhaustive
      N = _neighbours_definite_ex(res[i], p; callback, inv_dict, _invariants, use_mass, missing_mass, save_partial, save_path)
    elseif alg_type == :orbit
      N = _neighbours_definite_orbit(res[i], p; callback, inv_dict, _invariants, use_mass, missing_mass, save_partial, save_path)
    else
      N = _neighbours_definite_rand(res[i], p; rand_neigh, callback, inv_dict, _invariants, use_mass, missing_mass, save_partial, save_path)
    end
    if !is_empty(N)
      for M in N
        push!(tbv, true)
        push!(res, M)
      end
      use_mass && is_zero(missing_mass[]) && break
      if use_mass
        @v_do :ZZLatWithIsom perc = Float64(missing_mass[]//_mass) * 100
        @vprintln :ZZLatWithIsom 1 "Lattices: $(length(res)), Target mass: $(_mass). missing: $(missing_mass[]) ($(perc)%)"
      else
        @vprintln :ZZLatWithIsom "Lattices: $(length(res))"
      end
    end
  end
  return res, use_mass ? missing_mass[] : zero(QQ)
end

function enumerate_definite_genus(
    G::ZZGenus,
    alg_type::Symbol = :exhaustive;
    rand_neigh::Union{Int, Nothing} = nothing,
    invariant_func::Function = default_func,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true
  )
  _L = representative(G)
  neg = is_negative_definite(_L)
  if neg
    L = rescale(_L, -1)
  else
    L = _L
  end
  if use_mass
    s = automorphism_group_order(L)
    is_one(s*mass(G)) && return ZZLat[L], QQ(0)
    _missing_mass = mass(G) - 1//s
  end
  edg, mm = enumerate_definite_genus(ZZLat[L], alg_type; rand_neigh,
                                                         invariant_func,
                                                         save_partial,
                                                         save_path,
                                                         use_mass,
                                                         _missing_mass)

  while alg_type != :exhaustive && use_mass && !iszero(mm)
    edg, mm = enumerate_definite_genus(edg, alg_type; rand_neigh, invariant_func, save_partial, save_path, use_mass, _missing_mass = mm)
  end
  if neg
    map!(LL -> rescale(LL, -1), edg, edg)
  end
  return edg, mm
end

function enumerate_definite_genus(
    _L::ZZLat,
    alg_type::Symbol = :exhaustive;
    rand_neigh::Union{Int, Nothing} = nothing,
    invariant_func::Function = default_func,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true
  )
  neg = is_negative_definite(_L)
  if neg
    L = rescale(_L, -1)
  else
    L = _L
  end
  if use_mass
    s = automorphism_group_order(L)
    is_one(s*mass(L)) && return ZZLat[L], QQ(0)
    _missing_mass = mass(L) - 1//s
  end
  edg, mm = enumerate_definite_genus(ZZLat[L], alg_type; rand_neigh,
                                                         invariant_func,
                                                         save_partial,
                                                         save_path,
                                                         use_mass,
                                                         _missing_mass)

  while alg_type != :exhaustive && use_mass && !iszero(mm)
    edg, mm = enumerate_definite_genus(edg, alg_type; rand_neigh, invariant_func, save_partial, save_path, use_mass, _missing_mass = mm)
  end
  if neg
    map!(LL -> rescale(LL, -1), edg, edg)
  end
  return edg, mm
end
