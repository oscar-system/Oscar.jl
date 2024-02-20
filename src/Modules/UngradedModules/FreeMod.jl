###############################################################################
# FreeMod constructors
###############################################################################

@doc raw"""
    FreeMod(R::Ring, n::Int, name::VarName = :e; cached::Bool = false)

Construct a free module over the ring `R` with rank `n`.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod(R::Ring, n::Int, name::VarName = :e; cached::Bool = false) # TODO cached?
  return FreeMod{elem_type(R)}(n, R, [Symbol("$name[$i]") for i=1:n])
end

function FreeMod(R::Ring, names::Vector{String}; cached::Bool=false)
  return FreeMod{elem_type(R)}(length(names), R, Symbol.(names))
end

function FreeMod(R::Ring, names::Vector{Symbol}; cached::Bool=false)
  return FreeMod{elem_type(R)}(length(names), R, names)
end

@doc raw"""
    free_module(R::MPolyRing, p::Int, name::VarName = :e; cached::Bool = false)
    free_module(R::MPolyQuoRing, p::Int, name::VarName = :e; cached::Bool = false)
    free_module(R::MPolyLocRing, p::Int, name::VarName = :e; cached::Bool = false)
    free_module(R::MPolyQuoLocRing, p::Int, name::VarName = :e; cached::Bool = false)

Return the free $R$-module $R^p$, created with its basis of standard unit vectors.

The string `name` specifies how the basis vectors are printed. 

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> FR = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> x*FR[1]
x*e[1]

julia> P = ideal(R, [x, y, z]);

julia> U = complement_of_prime_ideal(P);

julia> RL, _ = localization(R, U);

julia> FRL = free_module(RL, 2, "f")
Free module of rank 2 over Localization of multivariate polynomial ring in 3 variables over QQ at complement of prime ideal (x, y, z)

julia> RL(x)*FRL[1]
x*f[1]

julia> RQ, _ = quo(R, ideal(R, [2*x^2-y^3, 2*x^2-y^5]));

julia> FRQ =  free_module(RQ, 2, "g")
Free module of rank 2 over RQ

julia> RQ(x)*FRQ[1]
x*g[1]

julia> RQL, _ = localization(RQ, U);

julia> FRQL =  free_module(RQL, 2, "h")
Free module of rank 2 over Localization of RQ at complement of prime ideal

julia> RQL(x)*FRQL[1]
x*h[1]
```
"""
free_module(R::MPolyRing, p::Int, name::VarName = :e; cached::Bool = false) = FreeMod(R, p, name, cached = cached)
free_module(R::MPolyQuoRing, p::Int, name::VarName = :e; cached::Bool = false) = FreeMod(R, p, name, cached = cached)
free_module(R::MPolyLocRing, p::Int, name::VarName = :e; cached::Bool = false) = FreeMod(R, p, name, cached = cached)
free_module(R::MPolyQuoLocRing, p::Int, name::VarName = :e; cached::Bool = false) = FreeMod(R, p, name, cached = cached)

#=XXX this cannot be as it is inherently ambiguous
  - FreeModule(R, n)
  - direct sum of rings, ie. a ring
  - set of n-th powers of R
thus the "category" needs to be set explicitly

=#


function (F::FreeMod)()
  return FreeModElem(sparse_row(base_ring(F)), F)
end

# by the magic of @show_name, this function will ensure that a 
# un-named free module over a named ring X will acquire the name 
# X^r
function AbstractAlgebra.extra_name(F::FreeMod)
  s = AbstractAlgebra.get_name(base_ring(F))
  if s !== nothing
    return "$s^$(rank(F))"
  end
  return nothing
end

function show(io::IO, F::FreeMod)
  @show_name(io, F)
  @show_special(io, F)
  compact = get(io, :compact, false)
  io_compact = IOContext(io, :compact => true)
  if is_graded(F)
      if !compact
        print(io, "Graded free module ")
      end
      i = 1
      while i <= dim(F)
          d = F.d[i]
          j = 1
          while i+j <= dim(F) && d == F.d[i+j]
              j += 1
          end
          print(io_compact, base_ring(F), "^$j")
          print(io_compact, "(", -d, ")")
          if i+j <= dim(F)
              print(io, " + ")
          end
          i += j
      end

      if rank(F)==0
        print(io_compact, base_ring(F), "^0")
      end

      if !compact
          print(io," of rank $(rank(F)) over ")
          print(io_compact, base_ring(F))
      end
  else
      if !compact
          #Todo: Use once the printing of rings is fixed
          #print(io_compact, "Free module ", base_ring(F), "^$(F.n) of rank $(F.n) over ")
          print(io_compact, "Free module of rank $(F.n) over ")
          print(io_compact, F.R)
      else
          print(io_compact, base_ring(F), "^$(F.n)")
      end
  end
end

@doc raw"""
    rank(F::FreeMod)
    number_of_generators(F::AbstractFreeMod)
    dim(F::AbstractFreeMod)

Return the rank of `F`.
"""
dim(F::AbstractFreeMod) = rank(F)
rank(F::FreeMod) = F.n
number_of_generators(F::AbstractFreeMod) = rank(F)

@doc raw"""
    ==(F::FreeMod, G::FreeMod)

Return  `true` if `F` and `G` are equal, `false` otherwise.

Here, `F` and `G` are equal iff either 
- both modules are ungraded and their base rings, ranks, and names for printing the basis elements are equal, 
or else 
- both modules are graded, the above holds, and for each $i$, the degrees of the $i$-th basis elements are equal.
"""
function (==)(F::FreeMod, G::FreeMod)
  # two free modules are equal if the rank and the ring are
  # TODO it this enough or e.g. stored morphisms also be considered?
  is_graded(F) == is_graded(G) || return false
  if is_graded(F) && is_graded(G) 
    return F.R == G.R && F.d == G.d && F.S == G.S
  end
  return F.R == G.R && rank(F) == rank(G) && F.S == G.S
end

function hash(F::FreeMod, h::UInt)
  b = is_graded(F) ? (0x2d55d561d3f7e215 % UInt) : (0x62ca4181ff3a12f4 % UInt)
  h = hash(base_ring(F), h)
  h = hash(rank(F), h)
  h = hash(F.S, h)
  is_graded(F) && (h = hash(F.d, h))
  return xor(h, b)
end

@doc raw"""
    is_isomorphic(F::FreeMod, G::FreeMod)

Return  `true` if `F` and `G` are isomorphic as (graded) modules, `false` otherwise.

That is, either 
- both modules are ungraded and their base rings and ranks are equal, 
or else 
- both modules are graded, the above holds, and the multisets of the degrees of the basis elements are equal.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, [1,1,3,2]);

julia> G1 = graded_free_module(Rg, [1,1,2,3]);

julia> is_isomorphic(F, G1)
true

julia> G2 = graded_free_module(Rg, [1,1,5,6]);

julia> is_isomorphic(F, G2)
false
```
"""
function is_isomorphic(F::FreeMod, G::FreeMod)
  is_graded(F) == is_graded(G) || return false
  is_graded(F) && return MSet(F.d) == MSet(G.d)
  return F.R == G.R && rank(F) == rank(G)
end

@doc raw"""
    is_zero(F::AbstractFreeMod)

Return `true` if `F` is the zero module, `false` otherwise.
"""
function is_zero(F::AbstractFreeMod)
  return rank(F) == 0
end

@doc raw"""
    canonical_isomorphism(F::FreeMod{T}, G::FreeMod{T})

For `F` and `G` which have equal rank (otherwise an error is thrown)
return the canonical isomorphism, that is map the i-th basis vector of the
canonical basis of `F` to the i-th basis vector of the canonical basis of `G`.
"""
function canonical_isomorphism(F::FreeMod{T}, G::FreeMod{T}) where T
  if F == G
     return hom(F, G, [G[i] for i in 1:ngens(G)])
  end
  @assert is_isomorphic(F, G)
  if is_graded(F) && is_graded(G)
    b = get_multiset_bijection(F.d, G.d, true)
    if length(b)==1
      b1 = b[1]
      return hom(F, G, [G[b1[i]] for i in 1:length(b1)])
    else
      error("there is no canonical isomorphism")
    end
    h = hom(F, G, [G[b[i]] for i in 1:length(b)])
    return h
  end
  if is_graded(F) != is_graded(G)
    error("there is no canonical isomorphism")
  end
  return hom(F, G, gens(G))
end

function isomorphism(F::FreeMod{T}, G::FreeMod{T}) where T
  @assert is_isomorphic(F, G)
  (!is_graded(F) && !is_graded(G)) && (return hom(F, G, gens(G)))
  if is_graded(F) != is_graded(G)
    error("there is no isomorphism")
  end
  b = get_multiset_bijection(F.d, G.d)
  h = hom(F, G, [G[b[i]] for i in 1:length(b)])
  return h
end

@doc raw"""
    ambient_free_module(F::FreeMod)

Just return `F`. This function exists only for compatibility reasons.
"""
function ambient_free_module(F::FreeMod)
  return F
end

@doc raw"""
    ambient_module(F::FreeMod, task = :none)

Just return `F`. This function exists only for compatibility reasons.

If task is set to :with_morphism, return also the identity map.
"""
function ambient_module(F::FreeMod, task = :none)
  if task == :none
    return F
  else
    return F, identity_map(F)
  end
end


@doc raw"""
    ambient_representatives_generators(F::FreeMod)

Return the generators of `F`. This function exists only for 
compatibility reasons.
"""
function ambient_representatives_generators(F::FreeMod)
  return gens(F)
end

rels(F::FreeMod) = elem_type(F)[]

