###############################################################################
#
#  Adaptation of the Hecke code to call Magma `LineOrbits`.
#  
#  Add also the possibility to run over the database of genera.
#
###############################################################################

##### Magma line orbits

#function magma_line_orbits(g::Vector{T}) where T <: MatElem
#  d = nrows(g[1])
#  F = base_ring(g[1])
#  V = vector_space(F, d)
#  p, e = characteristic(F), degree(F)
#  str = "K := GF($p, $e); G := MatrixGroup<$d, K | "
#  for m in g
#    mm = "[" * split(string([m[i,j] for i in 1:nrows(m) for j in 1:ncols(m)]), '[')[2]
#    str *= mm
#    str *= ", "
#  end
#  str = str[1:end-2]*" >; O := OrbitsOfSpaces(G, 1); Sprint([[[M[k] : k in [1..NumberOfColumns(M)]] : M in Basis(L[2])] : L in O])"
#  o = MagmaCall.interact() do stdout
#    MagmaCall.putcmd(stdout, str)
#    MagmaCall.readtotoken(String, stdout, missing) |> Meta.parse |> eval
#  end
#  orb = Vector{elem_type(F)}[F.(bas[1]) for bas in o]
#  return orb
#end

###### Neihbours

function _neighbours(
  L::ZZLat,
  p::ZZRingElem,
  algorithm::Symbol = :orbit;
  rand_neigh::Int=10,
  callback::Function=(M -> M != L),
  inv_dict=Ref(Dict()),
  _invariants::Function=(M -> []),
  save_partial::Bool=false,
  save_path::Union{Nothing, String}=nothing,
  use_mass::Bool=true,
  missing_mass::Base.RefValue{QQFieldElem}=Ref{QQFieldElem}(-1),
  vain::Base.RefValue{Int}=Ref{Int}(0),
  stop_after::IntExt=inf,
  max::IntExt=inf
)
  @assert !save_partial || !isnothing(save_path)
  bad = is_divisible_by(numerator(det(L)), p)
  even = is_even(L)
  K = GF(p)
  @assert algorithm in [:orbit, :random, :spinor]

  flag = (p == 2 && !even)
  if flag
    L0 = L
    L0 = even_sublattice(L)
    L0toL = solve(basis_matrix(L), basis_matrix(L0))
  else
    L0 = L
  end
  form0 = gram_matrix(L0)

  if even
    mqf = 2*p
    m = 2*p^2
  else
    mqf = p == 2 ? 4 : p
    m = p^2
  end

  if use_mass
    __mass = missing_mass[]
  end

  if algorithm == :orbit
    _, gene = isometry_group_smart(L) # For using Magma
    if flag
      bL0 = basis_matrix(L0)
      ibL0 = inv(bL0)
      gensp = dense_matrix_type(K)[map_entries(K, bL0*g*ibL0) for g in gene]
    else
      bL = basis_matrix(L)
      ibL = inv(bL)
      gensp = dense_matrix_type(K)[map_entries(K, bL*g*ibL) for g in gene]
    end
    @hassert :ZZLatWithIsom 3 !isempty(gensp)
############################## MAGMA ##########################################
    # orbs = try
    #          line_orbits_magma(gensp)
    #        catch
    # 	       Vector{elem_type(K)}[orb[1] for orb in Hecke.line_orbits(gensp)]
    #        end
############################## Hecke ##########################################
    orbs = Vector{elem_type(K)}[orb[1] for orb in Hecke.line_orbits(gensp)]
###############################################################################
    maxlines = length(orbs)
    stop_after = inf
    @vprintln :ZZLatWithIsom 3 "$(maxlines) orbits of lines to try"
  else
    P = Hecke.enumerate_lines(K, rank(L))
    maxlines = algorithm == :random ? min(rand_neigh, length(P)) : length(P)
    @vprintln :ZZLatWithIsom 3 "Try $(maxlines) random lines"
  end

  result = typeof(L)[]

  for i in 1:maxlines
    vain[] > stop_after && break
    if algorithm == :orbit
      x = orbs[i]
    elseif algorithm == :random
      x = rand(P)
    else
      x = Hecke.next(P)
    end
    w0 = matrix(QQ, 1, rank(L0), ZZRingElem[lift(ZZ, k) for k in x])
    a = numerator(only(w0*form0*transpose(w0)))
    if !is_divisible_by(a, mqf)
      vain[] += 1
      continue
    end

    lifts = typeof(w0)[]
    if p == 2
      bo = is_divisible_by(a, 8)
      if even && !bo
        w = w0*basis_matrix(L0)
        if is_zero(mod(divisibility(L, w), p)) # L_{w, 2} == L iff w lies in 2*L^#
          vain[] += 1
          continue
        end
        Hecke.make_admissible!(w0, form0, m, K, a)
        push!(lifts, w0)
      elseif even && bo # `w0` is admissible so it is good
        push!(lifts, w0)
      elseif !even && bo # Another corner case: `wL` is admissible but if `L_{wL, 2}` is even then the neighbour is even, and we want an odd one
        wL = w0*L0toL
        if is_even(Hecke.prime_dual(L, wL, p))
          vain[] += 1
          continue
        end
        push!(lifts, wL)
      else
        w = w0*basis_matrix(L0)
        wL = w0*L0toL
        # Here `w` is admissible of square 4 mod 8 so w/2 has odd square. Hence
        # `L_{w, 2} == L0` is allowed.
        #
        # If `L0_{w, 2} != L0`, we could allow another vector congruent to `w`
        # modulo 2*L0 with square divisible by 8. It is going to be another
        # admissible vector, and the neighbour might not be isometric to the one
        # obtained using `w`. The existence of such a vector is ensured only if
        # there exists a vector in `L0` with odd product with `w`, i.e. if
        # `L0_{w, 2} != L0`
        if is_zero(mod(divisibility(L0, w), p)) # L0_{w, 2} == L0 iff w lies in 2*L0^#
          push!(lifts, wL)
        else
          push!(lifts, wL)
          Hecke.make_admissible!(w0, form0, ZZ(8), K, a)
          push!(lifts, w0*L0toL)
        end
      end
    else
      if !is_divisible_by(a, m)
        w = w0*basis_matrix(L0)
        bad && is_zero(mod(divisibility(L, w), p)) && continue # In this case we cannot make w admissible because w lies in p*L^#
        Hecke.make_admissible!(w0, form0, m, K, a)
      end
      push!(lifts, w0)
    end

    for v in lifts
      LL = lll(Hecke.neighbour(L, v, p))
      @hassert :ZZLatWithIsom 3 is_locally_isometric(LL, L, p) # Should always hold by the neighbour construction

      keep = callback(LL)
      if !keep
        vain[] += 1
        continue
      end

      vain[] = Int(0)
      @vprintln :ZZLatWithIsom 3 "Keep an isometry class"
      if algorithm != :spinor
        invLL = _invariants(LL)
	if haskey(inv_dict[], invLL)
	  push!(inv_dict[][invLL], LL)
        else
	  inv_dict[][invLL] = ZZLat[LL]
        end
        push!(result, LL)

        if save_partial
          save_lattice(LL, save_path)
        end

        if use_mass
          s = isometry_group_order(LL)
	  sub!(__mass, __mass, 1//s)
          is_zero(__mass) && return result
        end

        length(result) == max && return result
      else
        return ZZLat[LL] # For :spinor we just want one neighbour
      end
    end
  end
  return result
end

##### Unique representative per isometry class
function __unique_iso_class!(A::Vector{ZZLat})
  isempty(A) && return A
  idxs = eachindex(A)

  # We keep the first lattice, so we start iterating from the second state
  it = iterate(idxs, 1)

  # Count the number of elements eventually kept
  # The function will iteratively stack from the start all the new
  # isometry classes to keep
  count = 1
  for x in Iterators.drop(A, 1)
    if all(y -> !is_isometric_smart(x, y), A[1:count])
      # In that case, x represents a new isometry class, so we stack it
      # after the first isometry classes kept at the beginning of `A`
      A[it[1]] = x
      count += 1
      it = iterate(idxs, it[2])
    end
  end
  # Now the lattices to keep are at the first `count` entries of `A`
  resize!(A, count)::typeof(A)
end

###### Default invariant function

@attr function _default_invariant_function(L::ZZLat)
  m = minimum(L)
  rlr, _ = root_lattice_recognition(L)
  kn = kissing_number(L)::Int
  ago = isometry_group_order(L)::ZZRingElem
  return (m, rlr, kn, ago)
end

###### Enumerate definite genus

function __enumerate_definite_genus(
    known::Vector{ZZLat},
    algorithm::Symbol = :default;
    rand_neigh::Int=10,
    distinct::Bool=false,
    invariant_function::Function=_default_invariant_function,
    save_partial::Bool=false,
    save_path::Union{String, Nothing}=nothing,
    use_mass::Bool=true,
    _missing_mass::Union{QQFieldElem, Nothing}=nothing,
    vain::Base.RefValue{Int}=Ref{Int}(0),
    stop_after::IntExt=inf,
    max::IntExt=inf
  )
  @req !save_partial || !isnothing(save_path) "No path mentioned for saving partial results"
  @req !is_empty(known) "Should know at least one lattice in the genus"
  @req all(LL -> genus(LL) == genus(known[1]), known) "Known lattices must be in the same genus"
  if algorithm != :default
    @req algorithm == :orbit || algorithm == :random "Only :random and :orbit algorithms are currently implemented"
  end

  @req !is_finite(max) || max > 0 "max must be infinite or positive"

  res = copy(known)
  !distinct && __unique_iso_class!(res)

  L, itk = Iterators.peel(res)
  inv_lat = invariant_function(L)
  inv_dict = Ref(Dict{typeof(inv_lat), Vector{ZZLat}}(inv_lat => ZZLat[L]))
  for N in itk
    inv_lat = invariant_function(N)
    if haskey(inv_dict[], inv_lat)
      push!(inv_dict[][inv_lat], N)
    else
      inv_dict[][inv_lat] = ZZLat[N]
    end
  end

  function _invariants(M::ZZLat)
    for (I, Z) in keys(inv_dict[])
      M in Z && return I
    end
    return invariant_function(M)
  end

  callback = function(M::ZZLat)
    any(isequal(M), res) && return false
    invM = _invariants(M)
    !haskey(inv_dict[], invM) && return true
    keep = all(N -> !is_isometric_smart(N, M), inv_dict[][invM])
    return keep
  end

  if use_mass
    _mass = mass(L)
    if isnothing(_missing_mass)
      found = sum(1//isometry_group_order(M) for M in res; init=QQ(0))
      missing_mass = Ref{QQFieldElem}(_mass-found)
      println(missing_mass[])
    else
      @hassert :GenRep 3 _missing_mass <= _mass
      missing_mass = Ref{QQFieldElem}(_missing_mass)
    end
  else
    missing_mass = Ref{QQFieldElem}(0)
  end

  r = rank(L)
  p = Hecke.smallest_kneser_prime(L)

  if algorithm == :default
    # Seems to be a reasonable bound for now, less than 1,000,000 lines
    if ndigits(divexact(ZZ(p)^r-1, p-1)) < 7
      algorithm = :orbit
    else
      algorithm = :random
    end
  end

  i = Int(0)
  while !is_zero(missing_mass[])
    i += 1
    N = _neighbours(res[i], p, algorithm; rand_neigh, callback, inv_dict, _invariants, use_mass, missing_mass, save_partial, save_path, vain, stop_after, max)

    if !is_empty(N)
      for M in N
        push!(res, M)
        if length(res) >= max
          return res, missing_mass[]
        end
      end
      use_mass && is_zero(missing_mass[]) && break
      if use_mass
        @v_do :ZZLatWithIsom 1 perc = Float64(missing_mass[]//_mass) * 100
        @vprintln :ZZLatWithIsom 1 "Lattices: $(length(res)), Target mass: $(_mass). missing: $(missing_mass[]) ($(perc)%)"
      else
        @vprintln :ZZLatWithIsom 1 "Lattices: $(length(res))"
      end
    elseif algorithm == :random && vain[] > stop_after
      break
    end
    if i == length(res)
      i = Int(0)
    end
  end
  return res, missing_mass[]
end

function __enumerate_definite_genus(
    L::ZZLat,
    algorithm::Symbol = :default;
    rand_neigh::Int = 10,
    invariant_function::Function = _default_invariant_function,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true,
    stop_after::IntExt = inf,
    max::IntExt = inf,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false
  )
  @req !save_partial || !isnothing(save_path) "No path mentioned for saving partial results"

  edg = ZZLat[]

  # We first compute representatives for each spinor genus
  spinor_genera = Hecke.spinor_genera_in_genus(L)
  @vprintln :ZZLatWithIsom 1 "$(length(spinor_genera)) spinor genera to enumerate"
  mm = QQ(0)
  # We enumerate each spinor genus separately.
  for M in spinor_genera
    vain = Ref{Int}(0)
    # The mass of a genus is "evenly distributed" among its spinor genera.
    # Hence, by dividing the mass of the genus by the number of spinor genera,
    # we know what portion of the mass the lattices of one spinor genus recover,
    # and we get a termination condition for each spinor genus.
    if use_mass
      mm = mass(M)//length(spinor_genera)
      s = isometry_group_order(M)
      sub!(mm, mm, 1//s)
      if is_zero(mm)
        push!(edg, M)
        continue
      end
    else
      @req 0 < stop_after < inf "Need to provide a finite positive value for stop_after if the mass is not used. Otherwise the algorithm may eventually never stops"
      mm = QQ(0)
    end
    _edg, mm = __enumerate_definite_genus(ZZLat[M], algorithm; rand_neigh,
                                                               invariant_function,
                                                               save_partial,
                                                               save_path,
                                                               use_mass,
                                                               _missing_mass=mm,
                                                               vain,
                                                               stop_after,
                                                               max=max-length(edg))
    while vain[] <= stop_after && length(edg) + length(_edg) < max
      use_mass && is_zero(mm) && break
      _edg, mm = __enumerate_definite_genus(_edg, algorithm; distinct=true,
                                                             rand_neigh,
                                                             invariant_function,
                                                             save_partial,
                                                             save_path,
                                                             use_mass,
                                                             _missing_mass=mm,
                                                             vain,
                                                             stop_after,
                                                             max=max-length(edg)-length(_edg))
    end
    append!(edg, _edg)
    length(edg) >= max && return edg
  end
  return edg, mm
end

function __enumerate_definite_genus(
    G::ZZGenus,
    algorithm::Symbol = :default;
    rand_neigh::Int=10,
    invariant_function::Function=_default_invariant_function,
    save_partial::Bool=false,
    save_path::Union{IO, String, Nothing}=nothing,
    use_mass::Bool=true,
    stop_after::IntExt=inf,
    max::IntExt=inf,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false
  )
  L = representative(G)
  @show G
  max == 1 && return ZZLat[L]
  return __enumerate_definite_genus(L, algorithm; rand_neigh,
                                                  invariant_function,
                                                  save_partial,
                                                  save_path,
                                                  use_mass,
                                                  stop_after,
                                                  max,
					                                        genusDB,
                                                  root_test)
end

function _smart_representatives(
  G::ZZGenus;
  genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
  root_test::Bool=false
)
  if !is_definite(G) || rank(G) <= 2
    return Hecke.representatives(G)
  end
  if !isnothing(genusDB)
    haskey(genusDB, G) && return genusDB[G]
  end
  r = rank(G)
  if root_test
    bn = Float64[0.5, 0.28868, 0.1847, 0.13127, 0.09987, 0.08112, 0.06981, 0.06326,
	       0.06007, 0.05953, 0.06136, 0.06559, 0.07253, 0.08278, 0.09735, 0.11774,
	       0.14624, 0.18629, 0.24308, 0.32454, 0.44289, 0.61722, 0.87767, 1.27241]
    if r <= 24 && abs(det(G)) < inv(bn[Int(r)])^2
      return ZZLat[]
    end
  end
  mm, l = __enumerate_definite_genus(G; genusDB, root_test, stop_after=1000)
  if !iszero(mm)
    inv_lat = _default_invariant_function(l[1])
    inv_dict = Dict{typeof(inv_lat), Vector{ZZLat}}(inv_lat => ZZLat[l[1]])
    for N in edg[2:end]
      inv_lat = _default_invariant_function(N)
      if haskey(inv_dict, inv_lat)
        push!(inv_dict[inv_lat], N)
      else
        inv_dict[inv_lat] = ZZLat[N]
      end
    end
    q = next_prime(last(Hecke.primes_up_to(r+1)))
    Lf = integer_lattice_with_isometry(l[1])
    pos = is_positive_definite(Lf)
    # Looking for certain lattices with isometry
    while !iszero(mm)
      d = denominator(mm)
      if isone(d)
        p = last(Hecke.primes_up_to(q-1))
      else
        p = maximum(prime_divisors(d))
      end
      @show d, p
      q = p
      if p == 2
        interv = div(r, 2):-1:1
      else
        interv = reverse(p-1:p-1:r)
      end
      for k in interv
        if pos
          Ns = splitting_of_prime_power(Lf, Int(p), 1; eiglat_cond=Dict(1=>[r-k, r-k, 0], p=>[k, k, 0]), genusDB, root_test=false, check=false)
        else
          Ns = splitting_of_prime_power(Lf, Int(p), 1; eiglat_cond=Dict(1=>[r-k, 0, r-k], p=>[k, 0, k]), genusDB, root_test=false, check=false)
        end
        for Nf in Ns
          N = lll(lattice(Nf))
          invN = _default_invariant_function(N)
          if !haskey(inv_dict, invN)
            inv_dict[invN] = ZZLat[N]
	          push!(edg, N)
	          s = isometry_group_order(N)
	          sub!(mm, mm, 1//s)
          elseif all(M -> !is_isometric_smart(N, M), inv_dict[invN])
            push!(inv_dict[invN], N)
            push!(edg, N)
	          s = isometry_group_order(N)
	          sub!(mm, mm, 1//s)
          end
          is_zero(mm) && break
        end
        is_zero(mm) && break
      end
    end
  end
  if !isnothing(genusDB)
    gesnuDB[G] = l
  end
  return l
end

###############################################################################
#
#  Compact storage of a genus
#
###############################################################################

# M must be symmetric
function __get_half_gram(L::ZZLat)
  M = gram_matrix(L)
  str = "$(nrows(M))\n["
  for i in 1:nrows(M), j in i:ncols(M)
    str *= "$(M[i,j]),"
  end
  str = str[1:end-1]*"]"
  if isdefined(L, :automorphism_group_order)
    str *= "\n$(L.automorphism_group_order))"
  end
  return str
end

# V is a list of gram matrices of the lattices in the genus
# f is a path to a numbered dir in the correct rank folder
function _save_genus(V::Vector{ZZLat}, f::String)
  for i in 1:length(V)
    p = f*"/lat_$(i).txt"
    touch(p)
    _f = open(p, "w")
    Base.write(_f, __get_half_gram(V[i]))
    close(_f)
  end
  return nothing
end

function _add_to_db!(f::String, V::Vector{ZZLat})
  @assert !isempty(V)
  G = genus(V[1])
  r = rank(G)
  fr = f*"rank$r/"
  if !isdir(fr)
    mkdir(fr)
  end
  l = length(readdir(fr))
  n = 5-ndigits(l+1)
  frl = fr*"0"^n*"$(l+1)"
  mkdir(frl)
  _save_genus(V, frl)
  return nothing
end

##############################################################################
#
#  Isometry testing
#
###############################################################################

function isometry_group_order(L::ZZLat)
  # corner case
  if rank(L) == 0
    return one(ZZ)
  end

  if isdefined(L, :automorphism_group_order)
    return L.automorphism_group_order
  end
#################################### MAGMA ####################################
  # LL = is_negative_definite(L) ? rescale(L, -1) : L
  # s = MagmaCall.interact() do stdout
  #   MagmaCall.putcmd(stdout, _to_magma_lattice(LL, "L")*"; G := AutomorphismGroup(L); Sprint(#G)")
  #   MagmaCall.readtotoken(String, stdout, missing) |> Meta.parse |> eval
  # end
  # L.automorphism_group_order = ZZ(s)
#################################### HECKE ####################################
  s = Hecke.automorphism_group_order(L)
###############################################################################
  return ZZ(s)
end

function is_isometric_smart(L::ZZLat, M::ZZLat; um=true)
  _default_invariant_function(L) != _default_invariant_function(M) && return false
#################################### MAGMA ####################################
  # LL = is_negative_definite(L) ? rescale(L, -1) : L
  # MM = is_negative_definite(M) ? rescale(M, -1) : M
  # b = MagmaCall.interact() do stdout
  #   MagmaCall.putcmd(stdout, _to_magma_lattice(LL, "L")*"; "*_to_magma_lattice(MM, "M")*"; Sprint(IsIsometric(L, M))")
  #   MagmaCall.readtotoken(String, stdout, missing) |> Meta.parse |> eval
  # end
#################################### HECKE ####################################
  b = Hecke.is_isometric(L, M)
###############################################################################
  return b
end

#################################### MAGMA ####################################
# function _to_magma_lattice(L::ZZLat, name::String)
#   Lmat = gram_matrix(L)
#   mat = "[" * split(string([Lmat[i,j] for i in 1:nrows(Lmat) for j in 1:ncols(Lmat)]), '[')[2]
#   mat = replace(mat, "//" => "/")
#   str = "$name := LatticeWithGram(Matrix(Rationals(), $(rank(L)), $(rank(L)), $mat))"
#   return str
# end
################################################################################

function isometry_group_smart(L::ZZLat)
  # corner case
  if rank(L) == 0
    gene = QQMatrix[identity_matrix(QQ,degree(L))]
    return matrix_group(gene), gene
  end

  if !is_definite(L) && (rank(L) == 2)
    gene1 = automorphism_group_generators(L; ambient_representation = false)
    gene2 = automorphism_group_generators(L; ambient_representation = true)
    return matrix_group(gene2), gene1
  end

  @req is_definite(L) "Lattice must be definite or of rank at most 2"
#################################### MAGMA ####################################  
  # if !isdefined(L, :automorphism_group_generators)
  #   G = gram_matrix(L)
  #   LL = is_negative_definite(L) ? integer_lattice(; gram = -G) : integer_lattice(; gram = G)
  #   _gene = MagmaCall.interact() do stdout
  #     MagmaCall.putcmd(stdout, _to_magma_lattice(LL, "L")*"; G := AutomorphismGroup(L); Gene := Generators(G); Sprint([[[M[j,k] : k in [1..NumberOfColumns(M)]] : j in [1..NumberOfRows(M)]] : M in Gene])")
  #     MagmaCall.readtotoken(String, stdout, missing) |> Meta.parse |> eval
  #   end
  #   V = ambient_space(L)
  #   B = basis_matrix(L)
  #   B2 = orthogonal_complement(V, B)
  #   C = vcat(B, B2)
  #   gene1 = QQMatrix[matrix(QQ, length(g), length(g), reduce(vcat, g)) for g in _gene]
  #   gene = QQMatrix[inv(C)*block_diagonal_matrix([m, identity_matrix(QQ, nrows(B2))])*C for m in gene1]
  #   aut = matrix_group(gene)
  # else
  #   _gene = L.automorphism_group_generators
  #   V = ambient_space(L)
  #   B = basis_matrix(L)
  #   B2 = orthogonal_complement(V, B)
  #   C = vcat(B, B2)
  #   gene1 = map(m -> map_entries(QQ, m), _gene)
  #   gene = QQMatrix[inv(C)*block_diagonal_matrix([m, identity_matrix(QQ, nrows(B2))])*C for m in gene1]
  #   aut = matrix_group(gene)
  # end
#################################### Hecke ####################################
  _gene = automorphism_group_generators(L)
  V = ambient_space(L)
  B = basis_matrix(L)
  B2 = orthogonal_complement(V, B)
  C = vcat(B, B2)
  gene1 = map(m -> map_entries(QQ, m), _gene)
  gene = QQMatrix[inv(C)*block_diagonal_matrix([m, identity_matrix(QQ, nrows(B2))])*C for m in gene1]
  aut = matrix_group(gene)
###############################################################################
  return aut, gene
end


############################### MAGMA #########################################
# function orbit_representatives_and_stabilizers_magma(G::MatrixGroup{E}, k::Int) where E <: FinFieldElem
#   g = gens(G)
#   d = degree(G)
#   F = base_ring(G)
#   V = vector_space(F, d)
#   p, e = characteristic(F), degree(F)
#   str = "K := GF($p, $e); G := MatrixGroup<$d, K | "
#   for m in g
#     mm = "[" * split(string([m[i,j] for i in 1:nrows(m) for j in 1:ncols(m)]), '[')[2]
#     str *= mm
#     str *= ", "
#   end
#   str = str[1:end-2]*" >; O := OrbitsOfSpaces(G, $k); S1 := Sprint([[[[M[j,k] : k in [1..NumberOfColumns(M)]] : j in [1..NumberOfRows(M)]] : M in Generators(Stabilizer(G, L[2]       ))] : L in O]); S2 := Sprint([[[M[k] : k in [1..NumberOfColumns(M)]] : M in Basis(L[2])] : L in O]); [S1, S2]"
#   o = MagmaCall.interact() do stdout
#     MagmaCall.putcmd(stdout, str)
#     MagmaCall.readtotoken(String, stdout, missing) |> Meta.parse |> eval
#   end
#   L1, L2 = o
#   stabs = Vector{elem_type(G)}[elem_type(G)[G(matrix(F, d, d, reduce(vcat, v))) for v in bas] for bas in L1]
#   stabs = [sub(G, bas)[1] for bas in stabs]   
#   orb = Vector{elem_type(V)}[elem_type(V)[V(F.(v)) for v in bas] for bas in L2]
#   orb = [sub(V, bas)[1] for bas in orb]
#   return [(orb[i], stabs[i]) for i in 1:length(orb)]
# end    
###############################################################################   
