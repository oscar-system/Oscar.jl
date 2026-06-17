include("exports.jl")
include("types.jl")

# TODO make docstrings prettier; specify R and J matrices
@doc raw"""
    modular_subgroup(s::PermGroupElem, t::PermGroupElem)

This function constructs a ModularGroup object corresponding to the finite-index subgroup
of SL2(Z) described by the permutations s and t. This constructor tests if the given
permutations actually describe the coset action of the matrices S = [0 -1; 1 0], T = [1 1; 0 1]
by checking that they act transitively and satisfy the relations
s^4 = (s^3 * t)^3 = s^2 * t * s^-2 * t^-1 = 1.

# Examples
```jldoctest
julia> s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
(1,2)(3,4)(5,6)(7,8)(9,10)

julia> t = cperm([1,4], [2,5,9,10,8], [3,7,6])
(1,4)(2,5,9,10,8)(3,7,6)

julia> G = modular_subgroup(s, t)
Modular subgroup of index 10
```
"""
function modular_subgroup(s::PermGroupElem, t::PermGroupElem)
  if !defines_coset_action_s_t(s, t)
    throw(ArgumentError("s and t do not describe the action of the generators S and T on the cosets of a finite-index subgroup of SL(2,Z)"))
  end
  return ModularGroup(s, t, s^-1*t^-1*s, s^-1*t^-1)
end

function Base.:(==)(G::ModularGroup, H::ModularGroup)
  return index(G) == index(H) && issubset(G, H)
end

function Base.hash(G::ModularGroup, h::UInt)
  return hash(G.s, hash(G.t, h))
end

function Base.show(io::IO, G::ModularGroup)
  idx = index(G)
  print(io, "Modular subgroup of index $(idx)")
end

@doc raw"""
    s_right_action(G::ModularGroup)

Returns the permutation describing the action of the matrix S on the right cosets of G.
"""
function s_right_action(G::ModularGroup)
  return G.s
end

@doc raw"""
    t_right_action(G::ModularGroup)

Returns the permutation describing the action of the matrix T on the right cosets of G.
"""
function t_right_action(G::ModularGroup)
  return G.t
end

@doc raw"""
    r_right_action(G::ModularGroup)

Returns the permutation describing the action of the matrix T on the right cosets of G.
"""
function r_right_action(G::ModularGroup)
  return G.r
end

@doc raw"""
    j_right_action(G::ModularGroup)

Returns the permutation describing the action of the matrix T on the right cosets of G.
"""
function j_right_action(G::ModularGroup)
  return G.j
end

function defines_coset_action_s_t(s::PermGroupElem, t::PermGroupElem)
  if !(isone(s^4) && isone((s^3*t)^3) && isone(s^2*t*s^-2*t^-1))
    return false
  end
  idx = maximum([maximum(moved_points(s); init=1), maximum(moved_points(t); init=1)])
  return is_transitive(permutation_group(idx, [s, t]))
end

function index(G::ModularGroup)
  return maximum([maximum(moved_points(G.s); init=1), maximum(moved_points(G.t); init=1)])::Int
end

const _SL2Z_FP_CACHE = Ref{Any}(nothing)

# cache the SL2Z presentation so that it can be reused consistently
# (especially in testing)
function _SL2Z_fp()
  if _SL2Z_FP_CACHE[] === nothing
    F = free_group(["S", "T"])
    S, T = gens(F)

    SL2Z, _ = quo(F, [S^4, (S^3*T)^3, S^2*T*S^-2*T^-1])

    S, T = gens(SL2Z)

    _SL2Z_FP_CACHE[] = (SL2Z, S, T)
  end

  return _SL2Z_FP_CACHE[]::Tuple{FPGroup, FPGroupElem, FPGroupElem}
end

function word_gens(G::ModularGroup)

  SL2Z, _, _ = _SL2Z_fp()
  P, _ = sub(symmetric_group(index(G)), [G.s, G.t])

  phi = hom(SL2Z, P, [G.s, G.t])

  # TODO: might need to check whether this is efficient enough for large index
  Hperm, _ = stabilizer(P, 1)

  H, inc = preimage(phi, Hperm)

  return [inc(h) for h in gens(H)]
end

function gens(G::ModularGroup)
  w_gens = word_gens(G)

  MatS = matrix(ZZ, [0 -1; 1 0])
  MatT = matrix(ZZ, [1  1; 0 1])

  M = matrix_group([MatS, MatT])
  MS, MT = gens(M)

  SL2Z, _, _ = _SL2Z_fp()
  phi = hom(SL2Z, M, [MS, MT])

  return [matrix(phi(w)) for w in w_gens]
end

function s_t_decomposition(M::ZZMatrix)
  if det(M) != 1
    throw(ArgumentError("Matrix needs to be in SL(2, Z)"))
  end

  MatS = matrix(ZZ, [0 -1; 1 0])
  MatT = matrix(ZZ, [1  1; 0 1])

  SL2Z, S, T = _SL2Z_fp()
  
  decomp = one(SL2Z)

  while M[2, 1] != 0
    k = div(M[2, 2], M[2, 1])
    decomp = S^-1 * T^k * decomp
    M = M * MatT^-k * MatS
  end

  # now M[2, 1] = 0 and since det(M) == 1, we have M == +-T^r where r = M[1, 2]
  if M[1, 1] == 1
    decomp = T^M[1, 2] * decomp
  else
    decomp = S^2 * T^-M[1, 2] * decomp
  end

  return decomp
end

function coset_action_of(A::ZZMatrix, G::ModularGroup)
  if det(A) != 1
     throw(ArgumentError("Matrix needs to be in SL(2, Z)"))
  end
  w = s_t_decomposition(A)
  P, _ = sub(symmetric_group(index(G)), [G.s, G.t])
  phi = hom(parent(w), P, [G.s, G.t])
  return phi(w)
end

function coset_right_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)
end

function Base.in(A::ZZMatrix, G::ModularGroup)
    if det(A) != 1
      return false
    end
    return coset_action_of(A, G)(1) == 1
end

function is_word_elm_of(w::FPGroupElem, G::ModularGroup)
  P, _ = sub(symmetric_group(index(G)), [G.s, G.t])
  phi = hom(parent(w), P, [G.s, G.t])
  return phi(w)(1) == 1
end

function Base.issubset(H::ModularGroup, G::ModularGroup)
  return all(A -> A in G, gens(H))
end
