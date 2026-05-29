@doc raw"""
    circuits(M::Matroid, S:: GroundsetType)

Return the list of circuits of the matroid contained in `S`. 

# Example
```jldoctest
julia> M = uniform_matroid(3,7);

julia> S = [1,3,4,6,7];

julia> circuits(M,S)
5-element Vector{Vector{Int64}}:
 [1, 3, 4, 6]
 [1, 3, 4, 7]
 [1, 3, 6, 7]
 [1, 4, 6, 7]
 [3, 4, 6, 7]
```
"""
function circuits(M::Matroid, S::T) where T<:GroundsetType
  @req issubset(S, matroid_groundset(M)) "The restriction set has to be a subset of the matroid's ground set"
  return circuits(restriction(M,S))
end

function _circuit(M::Matroid, S::T) where T<:GroundsetType
  C = circuits(M, S)
  @req !isempty(C) "S does not contain a circuit of the matroid"
  return C[1] 
end

@doc raw"""
    cocircuits(M::Matroid, S::T) where T<:GroundsetType

Return the list of cocircuits of matroid contained in `S`.

# Example
```jldoctest
julia> M = uniform_matroid(3,7);

julia> S = [1,3,4,6,7];

julia> cocircuits(M,S)
1-element Vector{Vector{Int64}}:
 [1, 3, 4, 6, 7]
```
"""
function cocircuits(M::Matroid, S::T) where T<:GroundsetType
  MD = dual_matroid(M)
  return circuits(MD,S)
end

function _cocircuit(M::Matroid, S::T) where T<:GroundsetType
  @req !iszero(rank(M)) "The rank of the matroid should be smaller than the cardinality of the ground set"
  return cocircuits(M,S)[1] 
end

@doc raw"""
    tutte_group(M::Matroid; char::Int=0)

Computes the Tutte group of a matroid `M` over a ring with characteristic `char`.
It should be noted, that the `char` only matters if it is two.
For more details see [DW89](@cite).

# Example
```jldoctest
julia> T = tutte_group(fano_matroid());

julia> ngens(T)
29
```
"""
function tutte_group(M::Matroid; char::Int=0)
  n = length(matroid_groundset(M))
  n > 128 && return _tutte_group_set(M, char)
  T = n<=8 ? UInt8 : n<=16 ? UInt16 : n<=32 ? UInt32 : n<=64 ? UInt64 : UInt128
  return _tutte_group(M, T, char)
end

function _tutte_group(M::Matroid, ::Type{T}, char::Int) where T <: Unsigned
  gs2num  = M.gs2num
  ebit(e) = T(1) << (gs2num[e] - 1)
  to_bits(s) = foldl((b, e) -> b | ebit(e), s; init=zero(T))

  B_sorted = sort([to_bits(b) for b in bases(M)])
  lookup(key) = searchsortedfirst(B_sorted, key)

  n       = length(matroid_groundset(M))
  gs_full = ~zero(T) >> (8*sizeof(T) - n)

  all_circ_vec    = circuits(M)
  all_cocirc_vec  = circuits(dual_matroid(M))
  all_circ_bits   = [to_bits(c) for c in all_circ_vec]
  all_cocirc_bits = [to_bits(c) for c in all_cocirc_vec]

  nb = length(B_sorted)
  v0 = zeros(Int, nb + 1); v0[end] = char == 2 ? 1 : 2
  relations = [v0]
  v_buf = zeros(Int, nb + 1)

  for X in nonbases(M)
    Xb = to_bits(X)
    ci = findfirst(c -> (c & Xb) == c, all_circ_bits)
    ci === nothing && continue
    findnext(c -> (c & Xb) == c, all_circ_bits, ci+1) !== nothing && continue
    di = findfirst(c -> (c & ~Xb & gs_full) == c, all_cocirc_bits)
    di === nothing && continue

    C  = all_circ_vec[ci];   Cb = [ebit(x) for x in C]
    D  = all_cocirc_vec[di]; Db = [ebit(x) for x in D]
    e = C[1]; eb = Cb[1]
    f = D[1]; fb = Db[1]

    for gi in 2:length(C)
      g = C[gi]; gb = Cb[gi]
      Ib   = Xb & ~(eb | gb)
      i_ef = lookup(Ib | eb | fb)
      i_gf = lookup(Ib | gb | fb)
      for hi in 2:length(D)
        h = D[hi]; hb = Db[hi]
        fill!(v_buf, 0)
        v_buf[i_ef]             =  1
        v_buf[lookup(Ib|gb|hb)] =  1
        v_buf[lookup(Ib|eb|hb)] = -1
        v_buf[i_gf]             = -1
        v_buf[end] = (e<f)+(f<g)+(g<h)+(h<e)
        push!(relations, copy(v_buf))
      end
    end
  end

  return abelian_group(matrix(ZZ, relations))
end

function _tutte_group_set(M::Matroid, char::Int)
  B = bases(M)
  gs = matroid_groundset(M)
  idx = Dict{Set{eltype(gs)}, Int}(Set(k) => i for (i,k) in enumerate(B))
  v0 = zeros(Int, length(B)+1); v0[end] = char == 2 ? 1 : 2
  relations = [v0]
  all_circ_vec   = circuits(M)
  all_cocirc_vec = circuits(dual_matroid(M))
  rkM = rank(M)
  v_buf = zeros(Int, length(B)+1)
  for X in nonbases(M)
    rank(M,X) == rkM-1 || continue
    Xset = Set(X)
    Yset = Set(setdiff(gs, X))
    ci = findfirst(c -> issubset(c, Xset), all_circ_vec)
    ci === nothing && continue
    di = findfirst(c -> issubset(c, Yset), all_cocirc_vec)
    di === nothing && continue
    C = copy(all_circ_vec[ci]); D = copy(all_cocirc_vec[di])
    e = popfirst!(C); f = popfirst!(D)
    for g in C
      I    = setdiff(Xset, (e, g))
      i_ef = idx[union(I, (e, f))]
      i_gf = idx[union(I, (g, f))]
      for h in D
        fill!(v_buf, 0)
        v_buf[i_ef]                =  1
        v_buf[idx[union(I,(g,h))]] =  1
        v_buf[idx[union(I,(e,h))]] = -1
        v_buf[i_gf]                = -1
        v_buf[end] = (e<f)+(f<g)+(g<h)+(h<e)
        push!(relations, copy(v_buf))
      end
    end
  end
  return abelian_group(matrix(ZZ, relations))
end

@doc raw"""
    is_tutte_realizable(M::Matroid)

Returns whether the matroid satisfies the Tutte realizability condition.
If `false`, this implies that the matroid cannot be realized over a field with characteristic other than `2`.
If `true`, we don't have a conclusive answer on realizability, since the Tutte
group only yields a necessary (and no sufficient) criterion for realizability
of `M`; see Corollary 1 in Section 3 of [DW89](@cite).

# Example
```jldoctest
julia> is_tutte_realizable(uniform_matroid(2,4))
true

julia> is_tutte_realizable(fano_matroid())
false
```
"""
function is_tutte_realizable(M::Matroid)
  T = tutte_group(M)
  return is_tutte_realizable(T)
end

function is_tutte_realizable(G::FinGenAbGroup)
  n=ngens(G)
  return !is_one(G[n])
end
