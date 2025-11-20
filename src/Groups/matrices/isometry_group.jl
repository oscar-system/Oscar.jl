"""
    isometry_group(L::AbstractLat; depth::Int = -1, bacher_depth::Int = 0) -> MatrixGroup

Return the group of isometries of the lattice `L`.

The transformations are represented with respect to the ambient space of `L`.

Setting the parameters `depth` and `bacher_depth` to a positive value may improve
performance. If set to `-1` (default), the used value of `depth` is chosen
heuristically depending on the rank of `L`. By default, `bacher_depth` is set to `0`.
"""
@attr matrix_group_type(S) function isometry_group(L::Hecke.AbstractLat{S}; depth::Int=-1, bacher_depth::Int=0) where S
  gens = automorphism_group_generators(L; depth, bacher_depth)
  G = matrix_group(gens)
  return G
end

@doc raw"""
    isometry_group(
      L::ZZLat;
      algorithm::Symbol=:default,
      depth::Int=-1,
      bacher_depth::Int=0,
      _set_nice_mono::Bool=true,
      _howell::Bool = true,
     ) -> MatrixGroup

Given an integer lattice $L$ which is definite or of rank 2, return the
isometry group $O(L)$ of $L$.

One can choose which algorithm to use to compute $O(L)$. For now, we
support the following algorithms:
- `:direct`: compute generators of $O(L)$ using an algorithm of
  Plesken--Souvignier;
- `:decomposition`: compute iteratively $O(L)$ by decomposing $L$ into
  $O(L)$-stable sublattices;
- `:default`: compute $O(L)$ using a heuristic to choose the best algorithm 
  between `:direct` and `:decomposition`.

Setting the parameters `depth` and `bacher_depth` to a positive value may
improve performance. If set to `-1` (default), the used value of `depth`
is chosen heuristically depending on the rank of `L`. By default,
`bacher_depth` is set to `0`.
"""
@attr QQMatrixGroup function isometry_group(
  L::ZZLat;
  algorithm::Symbol=:default,
  depth::Int=-1,
  bacher_depth::Int=0,
  _set_nice_mono::Bool=true,
  _howell::Bool = true,
)
  # G is represented w.r.t the basis of L
  G = Hecke._assert_has_automorphisms_ZZLat(L; algorithm, depth, bacher_depth, _set_nice_mono, _howell)
  Gamb = extend_to_ambient_space(L, G; check=false)
  if _set_nice_mono && is_definite(L) && rank(L) > 2
    set_order(Gamb, order(G))
    _sv = Hecke._short_vector_generators(L)
    sv = [i*basis_matrix(L) for i in _sv]
    _set_nice_monomorphism!(Gamb, sv)
  end 
  return Gamb
end

# Algorithm selection in `isometry_group` and `is_isometric_with_isometry` 
function _direct_is_faster(L::ZZLat)
    Llll = lll(L)
    G = gram_matrix(Llll)
    diagG = abs.(diagonal(G))
    ma = maximum(diagG)
    mi = minimum(diagG)
    r = rank(L)
    b =(r < 4 && ma <100*mi) || (r < 5 && ma <50*mi)|| (r < 6 && ma <25*mi)|| (r < 7 && ma <12*mi)|| (r < 8 && ma <9*mi) || (r < 9 && ma < 6*mi) || (r < 12 && ma < 4*mi)|| (ma < 2*mi)
    return b
end 
  
# We overwrite the function in Hecke in order that 
# Hecke has access to the better algorithms in Oscar
function Hecke._assert_has_automorphisms_ZZLat(L::ZZLat; 
                                               algorithm=:default,
                                               _howell::Bool = true,
                                               depth::Int=-1, 
                                               bacher_depth::Int=0,
                                               redo::Bool=false,
                                               _set_nice_mono::Bool=true
)
  # look in the cache
  if !redo && isdefined(L, :automorphism_group_generators)
    _gens = L.automorphism_group_generators
    G = matrix_group(_gens)
    if _set_nice_mono
      sv = Hecke._short_vector_generators(L)
      _set_nice_monomorphism!(G, sv)
    end
    return G
  end

  # corner cases
  @req rank(L) <= 2 || is_definite(L) "Lattice must be definite or of rank at most 2"
  if rank(L) <= 2
    Hecke.__assert_has_automorphisms(L; depth, bacher_depth, redo)
    _gens = L.automorphism_group_generators
    return matrix_group(_gens)
  end
  
  if algorithm == :default
    # Select an algorithm
    if _direct_is_faster(L::ZZLat)
      algorithm = :direct
    else
      algorithm = :decomposition
    end
  end
  
  if algorithm == :direct
    Hecke.__assert_has_automorphisms(L; depth, bacher_depth, redo)
    _gens = L.automorphism_group_generators
    G = matrix_group(_gens)
    if _set_nice_mono
      sv = Hecke._short_vector_generators(L)
      _set_nice_monomorphism!(G, sv)
    end
    set_order(G, L.automorphism_group_order) # computed in Hecke
  elseif algorithm == :decomposition
    G, _ = _isometry_group_via_decomposition(L; depth, bacher_depth, _howell)
    # fill the cache
    L.automorphism_group_order = order(G)
    L.automorphism_group_generators = ZZMatrix[change_base_ring(ZZ,matrix(g)) for g in gens(G)]
  else
    error("Unknown algorithm: for the moment, we only support :direct, :decomposition and :default")
  end
  return G
end

"""
    _isometry_group_via_decomposition(
      L::ZZLat;
      depth::Int=-1,
      bacher_depth::Int=0,
      direct::Bool=true,
      _set_nice_mono::Bool=true,
    ) -> MatrixGroup, Vector{QQMatrix}

Compute the group of isometries of the definite lattice `L` using an orthogonal
decomposition.
"""
function _isometry_group_via_decomposition(
  L::ZZLat;
  depth::Int=-1,
  bacher_depth::Int=0,
  _howell::Bool=true,
  _set_nice_mono::Bool=true
)
  # TODO: adapt the decomposition approach for AbstractLat
  # in most examples `direct=true` seems to be faster by a factor of 7
  # but in some examples it is also slower ... up to a factor of 15
  
  L = lattice(rational_span(L))
  if gram_matrix(L)[1,1] < 0
    L = rescale(L, -1)
  end

  if !_howell
    # need an integral lattice to work with discriminant groups
    d = denominator(scale(L))
    if d > 1
      L = rescale(L, d)
    end
  end

  # construct the sublattice M1 of L generated by the shortest vectors
  V = ambient_space(L)
  # for simplicity we work with the ambient representation
  M1, M1primitive, sv1 = Hecke._shortest_vectors_sublattice(L;check=false)
  basisM1prim = ZZ.(basis_matrix(M1primitive))
  # basically doubles the memory usage of this function
  # a more elegant way could be to work with the corresponding projective representation
  append!(sv1, [-v for v in sv1]) # given in the coordinates of L
    
  # select algorithm
  if _direct_is_faster(M1primitive)
    # compute O(M1primitive) directly in Hecke
    @vprintln :Isometry 3 "Computing orthogonal group in Hecke"
    Hecke.__assert_has_automorphisms(M1primitive; depth, bacher_depth) # avoid an infinite recursion
    O1 = matrix_group(M1primitive.automorphism_group_generators)
  else
    # first compute O(M1) then the stabiliser of M1primitive
    @vprintln :Isometry 3 "Computing orthogonal group of shortest sublattice in Hecke"
    Hecke.__assert_has_automorphisms(M1; depth, bacher_depth) # avoid an infinite recursion
    OM1 = matrix_group(M1.automorphism_group_generators)
    if M1primitive == M1
      _O1 = OM1
    else
      # M1 is generated by its shortest vectors only up to finite index
      @vprint :Isometry 3 "Computing overlattice stabilizers \n"
      svM1 = shortest_vectors(M1)
      append!(svM1, [-i for i in svM1])
      _set_nice_monomorphism!(OM1, svM1)
      @vtime :Isometry 3 _O1, _ = _overlattice_stabilizer(OM1, M1, M1primitive)
    end
    # transform to the basis of M1prim
    T = coordinates(basis_matrix(M1primitive), M1)
    invT = inv(T)
    O1 = matrix_group([ZZ.(T*matrix(g)*invT) for g in gens(_O1)])
    @hassert :Isometry 2 is_isometry_group(M1primitive, O1, false)
  end

  if rank(M1) == rank(L)
    B = basisM1prim 
    Binv = inv(B) 
    O1 = matrix_group([Binv*matrix(i)*B for i in gens(O1)])
    if _set_nice_mono
      _set_nice_monomorphism!(O1, sv1)
    end
    @hassert :Isometry 2 is_isometry_group(L, O1, false)
    return O1, sv1
  end
 
  # decompose as a primitive extension: M1primitive + M2 \subseteq L
  M2 = orthogonal_submodule(L, M1)

  # recursion
  @vprintln :Isometry 3 "Recursion\n"
  O2, sv2 = _isometry_group_via_decomposition(M2; _howell, depth, bacher_depth, _set_nice_mono)

  # go to to the basis of L
  basisM2 = ZZ.(basis_matrix(M2))
  sv = append!(sv1, [i*basisM2 for i in sv2])
  BB = vcat(basisM1prim,basisM2)
  
  # In what follows we compute the stabilizer of L in O1 x O2
  if _howell
    # Cook up O1 x O2
    r1 = rank(M1primitive)
    r2 = rank(M2) 
    I1 = identity_matrix(ZZ,r1)
    I2 = identity_matrix(ZZ,r2)
    gens12 = [block_diagonal_matrix(ZZMatrix[matrix(i),I2]) for i in gens(O1)] 
    append!(gens12, [block_diagonal_matrix(ZZMatrix[I1,matrix(i)]) for i in gens(O2)]) 
    O1timesO2 = matrix_group(gens12)
    B = vcat(basis_matrix(M1primitive),basis_matrix(M2))
    SL = lattice_in_same_ambient_space(L, B)
    _sv = Hecke._short_vector_generators(SL)
    _set_nice_monomorphism!(O1timesO2, _sv) # results in a vast speedup
    @vtime :Isometry 3 (_S,_) = _overlattice_stabilizer(O1timesO2, SL, L)
    BBs = solve_init(BB)
    # transform _S to the basis of L 
    _gens = [solve(BBs,matrix(i)*BB;side=:right) for i in gens(_S)]
    S = matrix_group(_gens)
  else
    O1amb = extend_to_ambient_space(M1primitive, O1; check=false)
    O2amb = extend_to_ambient_space(M2 ,O2; check=false)
    _sv1 = [i*basis_matrix(M1primitive) for i in Hecke._short_vector_generators(M1primitive)]
    _set_nice_monomorphism!(O1amb, _sv1)
    _sv2 = [i*basis_matrix(M2) for i in Hecke._short_vector_generators(M2)]
    _set_nice_monomorphism!(O1amb, _sv2)
    gensS = stabilizer_in_diagonal_action(L, M1primitive, M2, O1amb, O2amb; check=false, is_finite_known=(true, true))
    _S = matrix_group(gensS)
    SQ = restrict_to_lattice(L, _S)
    S = change_base_ring(ZZ,SQ)
  end

  if _set_nice_mono
    _set_nice_monomorphism!(S, sv)
  end
  @hassert :Isometry 2 is_isometry_group(L, S, false)
  return S, sv
end

function on_lattices(L::ZZLat, g::MatrixGroupElem{QQFieldElem,QQMatrix})
  V = ambient_space(L)
  return lattice(V, basis_matrix(L) * matrix(g); check=false)
end

"""
    on_vector(x::Vector{QQFieldElem}, g::MatrixGroupElem{QQFieldElem,QQMatrix})

Return `x*g`.
"""
function on_vector(x::Vector{QQFieldElem}, g::MatrixGroupElem{QQFieldElem,QQMatrix})
  return x*matrix(g)
end

"""
    _set_nice_monomorphism!(G::GAPGroup, nice_hom::GAPGroupHomomorphism)

Internally this sets a `NiceMonomorphism` for the underlying gap group.
It is assumed that the corresponding action homomorphism is injective.
No input checks whatsoever are performed.
"""
#
function _set_nice_monomorphism!(G::GAPGroup, nice_hom::GAPGroupHomomorphism)
  @req domain(nice_hom) === G "the domain of nice_hom must be G"
  f = GapObj(nice_hom)
  g = GapObj(G)
  set_is_finite(G, true)
  GAP.Globals.SetIsInjective(f, true) # fixes an infinite recursion
  GAP.Globals.SetIsHandledByNiceMonomorphism(g, true)
  GAP.Globals.SetNiceMonomorphism(g, f)
  return nothing
end

function _set_nice_monomorphism!(G::MatrixGroup{<:RingElem,T}, short_vectors::Vector{T}) where T<: Union{ZZMatrix,QQMatrix}
  nice_hom = _nice_hom!(G, copy(short_vectors))
  _set_nice_monomorphism!(G, nice_hom)
end 

function _set_nice_monomorphism!(G::MatrixGroup{S,T}, short_vectors::Vector{Vector{S}}) where {S<:Union{ZZRingElem, QQFieldElem}, T <: Union{ZZMatrix,QQMatrix}}
  R = base_ring(G)
  _short_vectors = T[matrix(R, 1, degree(G), i) for i in short_vectors]
  nice_hom = _nice_hom!(G, _short_vectors)
  _set_nice_monomorphism!(G, nice_hom)
end 

function _set_nice_monomorphism!(G::MatrixGroup{QQFieldElem,QQMatrix}, short_vectors::Vector{QQMatrix})
  nice_hom = _nice_hom!(G, copy(short_vectors))
  _set_nice_monomorphism!(G, nice_hom)
end 
  
# return g as permutation on X
# assumes that X is sorted
function _as_perm(w, g::T, X::Vector{T}) where T <: Union{ZZMatrix, QQMatrix}
  n = length(X)
  per = Vector{Int}(undef, n)
  i = 0
  for x in X
    i +=1
    mul!(w, x, g)
    j = searchsortedfirst(X, w,lt=Hecke._isless)
    per[i] = j
  end 
  return per
end 

# faster for big X
function _as_perm(w, g::T, X::IndexedSet{T}) where T <: Union{ZZMatrix, QQMatrix}
  n = length(X)
  per = Vector{Int}(undef, n)
  i = 0
  for i in 1:n
    mul!(w, X[i], g)
    per[i] = X(w)
  end 
  return per
end 

# sorts the _short_vectors !
function _nice_hom!(G::MatrixGroup{S, T}, _short_vectors::Vector{T}) where {S<:Union{ZZRingElem, QQFieldElem}, T<:Union{ZZMatrix, QQMatrix}}
  sort!(_short_vectors, lt=Hecke._isless)
  w = similar(_short_vectors[1])
  n = length(_short_vectors)
  Sn = symmetric_group(n)
  act_func(g) = perm(Sn, _as_perm(w, matrix(g), _short_vectors))
  return hom(G, Sn, act_func)
end 

# stabilizer of L in G < O(L)
function _overlattice_stabilizer(G::MatrixGroup{ZZRingElem,ZZMatrix}, S::ZZLat, L::ZZLat)
  _BL = coordinates(basis_matrix(L),S)
  n = denominator(_BL) 
  if n == 1
    # trivial nothing to do
    return G, hom(G,G,gens(G);check=false)
  end
  if is_prime(n)
    # up to 20% less allocations and slightly faster
    p = n
    BL = ZZ.(n*_BL)
    R = fpField(UInt(p))
    BLmod = change_base_ring(R, BL)
    r = rref!(BLmod)
    BLmod = BLmod[1:r,:]
    stab = stabilizer(G, BLmod, on_rref)
  else 
    BL = ZZ.(n*_BL)
    R,iR = residue_ring(ZZ, Int(n))
    BLmod = change_base_ring(R, BL)
    howell_form!(BLmod)
    stab = stabilizer(G, BLmod, on_howell_form)
  end 
  return stab
end 

function on_howell_form(M::zzModMatrix, g::MatrixGroupElem{ZZRingElem,ZZMatrix})
  Mg = M*matrix(g)
  howell_form!(Mg)
  return Mg
end 

function on_rref(M::fpMatrix, g::MatrixGroupElem{ZZRingElem,ZZMatrix})
  Mg = M*matrix(g)
  rref!(Mg)
  return Mg
end 


automorphism_group(L::Hecke.AbstractLat; kwargs...) = isometry_group(L; kwargs...)

orthogonal_group(L::Hecke.ZZLat; kwargs...) = isometry_group(L; kwargs...)

orthogonal_group(L::Hecke.QuadLat; kwargs...) = isometry_group(L; kwargs...)

unitary_group(L::Hecke.HermLat; kwargs...) = isometry_group(L; kwargs...)

@doc raw"""
    stable_orthogonal_group(
      L::ZZLat;
      kwargs...,
    ) -> MatrixGroup, GAPGroupHomomorphism

Given an integer lattice $L$ which is definite or of rank 2, return the
subgroup $O^\#(L)$ of the orthogonal group of $L$ consisting of isometries
acting trivially on the discriminant group of $L$.

The function first computes the orthogonal group of ``L``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> A5 = root_lattice(:A, 5);

julia> H, _ = stable_orthogonal_group(A5);

julia> order(H)
720
```
"""
function stable_orthogonal_group(
    L::ZZLat;
    kwargs...,
  )
  OL = orthogonal_group(L; kwargs...)
  return stable_subgroup(L, OL; check=false)
end

@doc raw"""
    special_orthogonal_group(
      L::ZZLat;
      kwargs...,
    ) -> MatrixGroup, GAPGroupHomomorphism

Given an integer lattice $L$ which is definite or of rank 2, return the
subgroup $SO(L)$ of the orthogonal group of $L$ consisting of isometries
with determinant ``1``.

The function first computes the orthogonal group of ``L``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> D5 = root_lattice(:D, 5);

julia> H, _ = special_orthogonal_group(D5);

julia> order(H)
1920
```
"""
function special_orthogonal_group(
    L::ZZLat;
    kwargs...,
  )
  OL = orthogonal_group(L; kwargs...)
  return special_subgroup(L, OL; check=false)
end

# We do not export this one, it is just a shortcut
@doc raw"""
    _special_stable_orthogonal_group(
      L::ZZLat;
      kwargs...,
    ) -> MatrixGroup, GAPGroupHomomorphism

Given an integer lattice $L$ which is definite or of rank 2, return the
subgroup $SO^\#(L)$ of the orthogonal group of $L$ consisting of isometries
acting trivially on the discriminant group of $L$ and of determinant ``1``.

The function first computes the orthogonal group of ``L``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> A6 = root_lattice(:A, 6);

julia> H, _ = Oscar._special_stable_orthogonal_group(A6);

julia> describe(H)
"A7"
```
"""
function _special_stable_orthogonal_group(
    L::ZZLat;
    kwargs...,
  )
  OL = orthogonal_group(L; kwargs...)
  return Oscar._special_stable_subgroup(L, OL; check=false)
end

function Hecke._is_isometric_with_isometry_definite(L1::ZZLat, 
                                                    L2::ZZLat; 
                                                    kwargs...)
   return _is_isometric_with_isometry_definite_via_decomposition(L1, L2; kwargs...)
end

function _is_isometric_with_isometry_definite_via_decomposition(L1::ZZLat, 
                                                                L2::ZZLat; 
                                                                depth::Int = -1, 
                                                                bacher_depth::Int = 0)
  L1 = lattice(rational_span(L1))
  L2 = lattice(rational_span(L2))
  
  # algorithm selection
  
  if _direct_is_faster(L1)
    b,fL = Hecke.__is_isometric_with_isometry_definite(L1, L2; depth, bacher_depth)
    @hassert :Isometry 3 fL*gram_matrix(L2)*transpose(fL) == gram_matrix(L1)
    return b, fL
  end
  
  # todo: if degree > rank go to lattice(rational_span) and transform back?
  rank(L1) != rank(L2) && return false, zero_matrix(QQ, 0, 0)
  # TODO: compute short vectors only once
  # TODO: start with the non-primitive one
  M1s, M1,_ = Hecke._shortest_vectors_sublattice(L1; check=false)
  M2s, M2,_ = Hecke._shortest_vectors_sublattice(L2; check=false)
  
  # avoid that Hecke enters a possibly expensive computation
  # the minimum is known anyways
  minimum(L1) == minimum(L2) || return false, zero_matrix(QQ, 0, 0)
  # algorithm selection
  if _direct_is_faster(M1)
    b, fM = Hecke.__is_isometric_with_isometry_definite(M1, M2;  depth, bacher_depth)
    b || return false, zero_matrix(QQ, 0, 0)
  else 
    b, fMs = Hecke.__is_isometric_with_isometry_definite(M1s, M2s;  depth, bacher_depth)
    b || return false, zero_matrix(QQ, 0, 0)
    OM1s,_ = _isometry_group_via_decomposition(M1s; depth, bacher_depth)
    # make it work for now. Go through hecke later?
    @vprint :Isometry 2 "computing orbit of an overlattice..."
    DM1s = discriminant_group(M1s)
    dM1s = discriminant_representation(M1s, OM1s; check=false, full=false, ambient_representation=false)
    OM1sq = image(dM1s)[1]
    to_gapM1s = get_attribute(OM1sq,:to_gap)
    B1 = basis_matrix(M1)
    J1gap = sub(codomain(to_gapM1s), [to_gapM1s(DM1s(B1[i,:])) for i in 1:nrows(B1)])[1]
    fMsinvB = inv(fMs)*basis_matrix(M1s)
    B2 = basis_matrix(M2)
    J2gap = sub(codomain(to_gapM1s), [to_gapM1s(DM1s(coordinates(B2[i,:],M2s)*fMsinvB)) for i in 1:nrows(B2)])[1]
    XM = orbit(OM1sq, on_subgroups, J1gap)
    b, tmp = is_conjugate_with_data(XM, J1gap, J2gap)
    b || return false, zero_matrix(QQ, 0, 0)
    @vprintln :Isometry 2 "done"
    # transform back to the bases of M1 and M2
    T1 = solve(basis_matrix(M1), basis_matrix(M1s); side=:left)
    T2 = solve(basis_matrix(M2), basis_matrix(M2s); side=:left)
    fM = inv(T1)*matrix(dM1s\tmp)*fMs*T2
  end
  @hassert :Isometry 3 fM*gram_matrix(M2)*transpose(fM) == gram_matrix(M1)
  
  # base case of the recursion 
  if rank(M1) == rank(L1)
    @assert degree(M1)==rank(M1)
    @assert degree(M2)==rank(M2)
    f = inv(basis_matrix(M1))*fM*basis_matrix(M2)
    @hassert :Isometry 1 lattice_in_same_ambient_space(L2, basis_matrix(L1)*f) == L2
    @hassert :Isometry 1 gram_matrix(ambient_space(L2), f) == gram_matrix(ambient_space(L1))
    return b, f
  end
  # recurse 
  N1 = orthogonal_submodule(L1, M1)
  N2 = orthogonal_submodule(L2, M2)
  b, fN = is_isometric_with_isometry(N1, N2;  depth, bacher_depth)
  b || return false, zero_matrix(QQ, 0, 0)
  @hassert :Isometry 3 fN*gram_matrix(N2)*transpose(fN) == gram_matrix(N1)
  #b, fN = _is_isometric_via_decomposition(N1, N2)
  

  
  # modify (fM,fN) so that it extends to L
  @vtime :Isometry 4 OM1,_ = _isometry_group_via_decomposition(M1; depth, bacher_depth)
  @vtime :Isometry 4 ON1,_ = _isometry_group_via_decomposition(N1; depth, bacher_depth)
  
  DM1 = discriminant_group(M1)
  dM1 = discriminant_representation(M1, OM1; check=false, full=false, ambient_representation=false)
  OM1q = image(dM1)[1]
  DN1 = discriminant_group(N1)
  dN1 = discriminant_representation(N1, ON1; check=false, full=false, ambient_representation=false)
  ON1q = image(dN1)[1]
  
  phi1, iHM1, iHN1 = glue_map(L1, M1,N1;_snf=true, check=false)
  phi2, iHM2, iHN2 = glue_map(L2, M2,N2;_snf=true, check=false)
  
  # TODO: early out if direct sum. 
  
  # Modify (fM,fN) so that fM(HM1) = HM2 and fN(HN1) = HN2
  HM1 = domain(phi1)
  HN1 = codomain(phi1)
  HM2 = domain(phi2)
  HN2 = codomain(phi2)
  
  to_gapM1 = get_attribute(OM1q,:to_gap)
  #to_oscarM = get_attribute(OM1q,:to_oscar)
  HM1gap = sub(codomain(to_gapM1), [to_gapM1(iHM1(i)) for i in gens(HM1)])[1]
  fMinvB = inv(fM)*basis_matrix(M1)
  HM2gap = sub(codomain(to_gapM1), [to_gapM1(DM1(coordinates(lift(i),M2)*fMinvB)) for i in gens(HM2)])[1]
  @vtime :Isometry 4 XM = orbit(OM1q, on_subgroups, HM1gap)
  @vtime :Isometry 4 (b, tmp) = is_conjugate_with_data(XM, HM1gap, HM2gap)
  b || return false, zero_matrix(QQ, 0, 0)
  
  fM = matrix(dM1\tmp)*fM
  fMB = fM*basis_matrix(M2)
  # confirm result
  @hassert :Isometry 3 all(coordinates(lift(i),M1)*fMB in cover(HM2) for i in gens(HM1))
  @hassert :Isometry 3 fM*gram_matrix(M2)*transpose(fM) == gram_matrix(M1)

  to_gapN1 = get_attribute(ON1q,:to_gap)
  #to_oscarN = get_attribute(ON1q,:to_oscar)
  HN1gap = sub(codomain(to_gapN1), [to_gapN1(iHN1(i)) for i in gens(HN1)])[1]
  fNinvB = inv(fN)*basis_matrix(N1)
  HN2gap = sub(codomain(to_gapN1), [to_gapN1(DN1(coordinates(lift(i),N2)*fNinvB)) for i in gens(HN2)])[1]
  XN = orbit(ON1q, on_subgroups, HN1gap)
  
  b, tmp = is_conjugate_with_data(XN, HN1gap, HN2gap)
  b || return false, zero_matrix(QQ, 0, 0)
  # fN is represented w.r.t. the bases of N1 and N2
  fN = matrix(dN1\tmp)*fN
  @hassert :Isometry 3 fN*gram_matrix(N2)*transpose(fN) == gram_matrix(N1)
  fNB = fN*basis_matrix(N2)
  # confirm result
  @hassert :Isometry 3 all(coordinates(lift(i),N1)*fNB in cover(HN2) for i in gens(HN1))
  
  fHM1_HM2 = hom(HM1, HM2, [HM2(coordinates(lift(i), M1)*fMB) for i in gens(HM1)])
  fHN1_HN2 = hom(HN1, HN2, [HN2(coordinates(lift(i), N1)*fNB) for i in gens(HN1)])
  # go in a square and obtain an isometry of HN1 
  # if the square commutes, g is trivial 
  # thus we check if we can make g trivial by modifying f
  _g = fHM1_HM2*phi2*inv(fHN1_HN2)*inv(phi1)
  OHM1 = orthogonal_group(HM1) 
  g = OHM1(_g)
  
  # setup the induced orthogonal groups
  stabHM1,inc_stabHM1 = stabilizer(OM1q, HM1gap, on_subgroups)
  stabHN1,inc_stabHN1 = stabilizer(ON1q, HN1gap, on_subgroups)
  stabHM1_on_HM1,res_stabM = restrict_automorphism_group(stabHM1, iHM1; check=false)
  stabHN1_on_HN1,res_stabN= restrict_automorphism_group(stabHN1, iHN1;check=false)
  phi_stabHN1_on_HN1_phiinv,_ = sub(OHM1, elem_type(OHM1)[OHM1(phi1*hom(i)*inv(phi1)) for i in gens(stabHN1_on_HN1)])
  
  
  # notation D:= UgV \subseteq G
  # find u,v with 1 = u g v
  G = OHM1
  U = stabHM1_on_HM1
  V = phi_stabHN1_on_HN1_phiinv
  
  # permutation groups are faster
  GtoGperm = isomorphism(PermGroup, G) # there are other ways to get a permutation group
  Gperm = codomain(GtoGperm)
  Uperm = GtoGperm(U)[1]
  Vperm = GtoGperm(V)[1]
  g1perm = one(Gperm)
  gperm = GtoGperm(g)
  D = double_coset(Uperm, gperm, Vperm)
  b, uperm, vperm = _decompose(D, g1perm)
  b || return false, zero_matrix(QQ, 0, 0)
  v = GtoGperm\vperm 
  u = GtoGperm\uperm
  @hassert :Isometry 3 u*v == g
  
  gMbar = inv(u) 
  gNbar = inv(v)
  @hassert :Isometry 3 one(G) == gMbar * g * gNbar 
  gNbar = stabHN1_on_HN1(inv(phi1)*hom(v)*phi1)
  f1 = inc_stabHM1\(res_stabM\gMbar)
  f2 = inc_stabHN1\(res_stabN\gNbar)
    
  fM = (dM1\f1)*fM
  fN = (dN1\f2)*fN
  @hassert :Isometry 3 fM*gram_matrix(M2)*transpose(fM) == gram_matrix(M1)
  @hassert :Isometry 3 fN*gram_matrix(N2)*transpose(fN) == gram_matrix(N1)
  # assemble the isometry 
  
  B1 = vcat(basis_matrix(M1),basis_matrix(N1))
  B2 = vcat(basis_matrix(M2),basis_matrix(N2))
  f = solve(B1, diagonal_matrix([fM,fN]); side=:right)*B2
  @hassert :Isometry 3 lattice_in_same_ambient_space(L2, basis_matrix(M1)*f) == M2
  @hassert :Isometry 3 lattice_in_same_ambient_space(L2, basis_matrix(N1)*f) == N2

  @hassert :Isometry 1 lattice_in_same_ambient_space(L2, basis_matrix(L1)*f) == L2
  @hassert :Isometry 1 gram_matrix(ambient_space(L2), f) == gram_matrix(ambient_space(L1))
  return b, f
end
