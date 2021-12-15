export SpecOpen, parent, gens, complement, patches, intersections, name, intersect, union, issubset, closure

export StructureSheafRing, variety, domain, OO

export StructureSheafElem, domain, representatives, patches, npatches, restrict

export SpecOpenMor, maps_on_patches, restriction

mutable struct SpecOpen{BRT, BRET, RT, RET, MST} <: Scheme{BRT, BRET}
  X::Spec{BRT, BRET, RT, RET, MST} # the ambient scheme
  gens::Vector{RET} # a list of functions defining the complement of the open subset

  # fields used for caching
  name::String
  patches::Vector{Spec{BRT, BRET, RT, RET, MST}}
  intersections::Dict{Tuple{Int}, Spec{BRT, BRET, RT, RET, MST}}
  complement::Spec{BRT, BRET, RT, RET, MST}

  function SpecOpen(X::Spec{BRT, BRET, RT, RET, MST}, f::Vector{RET}; name::String="") where {BRT, BRET, RT, RET, MST}
    for a in f
      parent(a) == base_ring(OO(X)) || error("element does not belong to the correct ring")
    end
    if length(name) > 0 
      U = new{BRT, BRET, RT, RET, MST}(OO(X), f)
      set_name!(U, name)
      return U
    end
    return new{BRT, BRET, RT, RET, MST}(X, f)
  end
end

parent(U::SpecOpen) = U.X
gens(U::SpecOpen) = U.gens

function Base.show(io::IO, U::SpecOpen)
  if isdefined(U, :name) 
    print(io, name(U))
    return
  end
  print(io, "complement of zero locus of $(gens(U)) in $(parent(U))")
end

function SpecOpen(
    X::Spec{BRT, BRET, RT, RET, MST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  base_ring(I) == localized_ring(OO(X)) || error("Ideal does not belong to the correct ring")
  f = [reduce(f, groebner_basis(localized_modulus(OO(X)))) for f in gens(I)]
  g = Vector{elem_type(base_ring(OO(X)))}()
  for a in f
    iszero(numerator(a)) || (push!(g, numerator(a)))
  end
  return SpecOpen(X, g)
end

function SpecOpen(
    X::Spec{BRT, BRET, RT, RET, MST},
    I::MPolyIdeal{RET}
  ) where {BRT, BRET, RT, RET, MST}
  return SpecOpen(X, localized_ring(OO(X))(I))
end

function complement(X::T, Z::T) where {T<:Spec}
  if !issubset(Z, X) 
    Z = intersect(X, Z)
  end
  return SpecOpen(X, modulus(OO(Z)))
end

SpecOpen(X::Spec) = SpecOpen(X, [one(base_ring(OO(X)))])

function complement(U::SpecOpen) 
  if !isdefined(U, :complement)
    I = radical(saturated_ideal(ideal(localized_ring(OO(parent(U))), gens(U))))
    U.complement = subscheme(parent(U), I)
  end
  return U.complement
end

function patches(U::SpecOpen)
  if !isdefined(U, :patches)
    X = parent(U)
    V = Vector{typeof(X)}()
    for f in gens(U)
      push!(V, hypersurface_complement(X, f))
    end
    U.patches = V
  end
  return U.patches
end

function intersections(U::SpecOpen)
  if !isdefined(U, :intersections)
    X = parent(U)
    V = patches(U)
    for i in 2:length(V)
      for j in 1:i-1
	U.intersections[(i,j)] = U.intersections[(j,i)] = intersect(V[i], V[j])
      end
    end
  end
  return U.intersections
end

function name(U::SpecOpen) 
  if isdefined(U, :name)
    return U.name
  end
  return "open subset of $(parent(U))"
end

function intersect(
    Y::Spec{BRT, BRET, RT, RET, MST}, 
    U::SpecOpen{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  X = parent(U)
  if !issubset(Y, X)
    Y = intersect(Y, X)
  end
  return SpecOpen(Y, gens(U))
end

function intersect(
    U::SpecOpen{BRT, BRET, RT, RET, MST},
    Y::Spec{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  return intersect(Y, U)
end

function union(U::T, V::T) where {T<:SpecOpen}
  parent(U) == parent(V) || error("the two open sets do not lay in the same ambient variety")
  return SpecOpen(parent(U), vcat(gens(U), gens(V)))
end

function issubset(
    Y::Spec{BRT, BRET, RT, RET, MST}, 
    U::SpecOpen{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  return one(OO(Y)) in ideal(OO(Y), gens(U))
end

function issubset(
    U::SpecOpen{BRT, BRET, RT, RET, MST},
    Y::Spec{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  for V in patches(U)
    issubset(V, Y) || return false
  end
  return true
end

function issubset(U::T, V::T) where {T<:SpecOpen}
  return issubset(complement(U), complement(intersect(V, parent(U))))
end

function ==(U::T, V::T) where {T<:SpecOpen}
  return issubset(U, V) && issubset(V, U)
end

function ==(
    U::SpecOpen{BRT, BRET, RT, RET, MST},
    Y::Spec{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  return issubset(U, Y) && issubset(Y, U)
end

function ==(
    Y::Spec{BRT, BRET, RT, RET, MST},
    U::SpecOpen{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  return U == Y
end

function closure(
    U::SpecOpen{BRT, BRET, RT, RET, MST},
    Y::Spec{BRT, BRET, RT, RET, MST} 
  ) where {BRT, BRET, RT, RET, MST}
  error("not implemented")
end


@Markdown.doc """
StructureSheafRing{BRT, BRET, RT, RET, MST}

The ring of regular functions 𝒪(X, U) on an open subset U of an 
affine scheme X.
"""
mutable struct StructureSheafRing{BRT, BRET, RT, RET, MST}
  variety::Spec{BRT, BRET, RT, RET, MST}
  domain::SpecOpen{BRT, BRET, RT, RET, MST}

  function StructureSheafRing(
      X::Spec{BRT, BRET, RT, RET, MST}, 
      U::SpecOpen{BRT, BRET, RT, RET, MST}
    ) where {BRT, BRET, RT, RET, MST}
    issubset(U, X) || error("open set does not lay in the variety")
    return new{BRT, BRET, RT, RET, MST}(X, U)
  end
end

variety(R::StructureSheafRing) = R.variety
domain(R::StructureSheafRing) = R.domain

OO(U::SpecOpen) = StructureSheafRing(parent(U), U)
OO(X::Spec{BRT, BRET, RT, RET, MST}, U::SpecOpen{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = StructureSheafRing(X, U)

function ==(R::T, S::T) where {T<:StructureSheafRing} 
  variety(R) == variety(S) || return false
  domain(S) == domain(R) || return false
  return true
end

elem_type(R::StructureSheafRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}= Type{StructureSheafElem{BRT, BRET, RT, RET, MST}}

@Markdown.doc """
StructureSheafElem{BRT, BRET, RT, RET, MST}

An element f ∈ 𝒪(X, U) of the ring of regular functions on 
an open set U of an affine variety X.
"""
mutable struct StructureSheafElem{BRT, BRET, RT, RET, MST}
  domain::SpecOpen{BRT, BRET, RT, RET, MST}
  representatives::Vector{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}

  function StructureSheafElem(
      U::SpecOpen{BRT, BRET, RT, RET, MST},
      f::Vector{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}
    ) where {BRT, BRET, RT, RET, MST}
    n = length(f)
    n == length(patches(U)) || error("the number of representatives does not coincide with the number of affine patches")
    for i in 1:n 
      OO(patches(U)[i])(lift(f[i])) # throws an error if conversion is not possible
    end
    return new{BRT, BRET, RT, RET, MST}(U, f)
  end
end

variety(f::StructureSheafElem) = parent(domain(f))
domain(f::StructureSheafElem) = f.domain
representatives(f::StructureSheafElem) = f.representatives
patches(f::StructureSheafElem) = patches(domain(f))
npatches(f::StructureSheafElem) = length(f.representatives)
getindex(f::StructureSheafElem, i::Int) = getindex(representatives(f), i)
parent(f::StructureSheafElem) = OO(variety(f), domain(f))

function restrict(
    f::StructureSheafElem{BRT, BRET, RT, RET, MST}, 
    V::Spec{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  for i in 1:length(representatives(f))
    if V == patches(domain(f))[i]
      return representatives(f)[i]
    end
  end
  issubset(V, domain(f)) || error("the set is not contained in the domain of definition of the function")
  if !is_open_embedding(V, parent(domain(f)))
    W = intersect(V, domain(f))
    g = StructureSheafElem(W, [OO(patches(W)[i])(f[i]) for i in 1:npatches(f)])
    return restrict(g, V)
  end
  error("not implemented")
end

@Markdown.doc """
maximal_extension(
  X::Spec{BRT, BRET, RT, RET, MST}, 
  f::AbstractAlgebra.Generic.Frac{RET}
) where {BRT, BRET, RT, RET, MST}

Returns the maximal extension of the restriction of f 
to a rational function on X on a maximal domain of 
definition U ⊂ X. 
"""
function maximal_extension(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {BRT, BRET, RT, RET, MST}
  a = numerator(f)
  b = denominator(f)
  W = localized_ring(OO(X))
  I = quotient(ideal(W, b) + localized_modulus(OO(X)), ideal(W, a))
  U = SpecOpen(X, I)
  g = [OO(V)(f) for V in patches(U)]
  return StructureSheafElem(U, g)
end

@Markdown.doc """
maximal_extension(
  X::Spec{BRT, BRET, RT, RET, MST}, 
  f::Vector{AbstractAlgebra.Generic.Frac{RET}}
) where {BRT, BRET, RT, RET, MST}

Returns the extension of the restriction of the fᵢ as a 
set of rational functions on X to a common maximal domain of 
definition U ⊂ X.
"""
function maximal_extension(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {BRT, BRET, RT, RET, MST}
  if length(f) == 0
    return SpecOpen(X), Vector{StructureSheafElem{BRT, BRET, RT, RET, MST}}()
  end
  a = numerator.(f)
  b = denominator.(f)
  W = localized_ring(OO(X))
  I = ideal(W, one(W))
  for p in f
    I = intersection(quotient(ideal(W, denominator(p)) + localized_modulus(OO(X)), ideal(W, numerator(p))), I)
  end
  U = SpecOpen(X, I)
  return U, [StructureSheafElem(U, [OO(V)(a) for V in patches(U)]) for a in f]
end

@Markdown.doc """
SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2}

Morphisms f : U → V of open sets U ⊂ X and V ⊂ Y of affine schemes.
These are stored as morphisms fᵢ: Uᵢ→ Y on the affine patches 
Uᵢ of U.
"""
mutable struct SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2}
  domain::SpecOpen{BRT, BRET, RT, RET, MST1}
  codomain::SpecOpen{BRT, BRET, RT, RET, MST2}
  maps_on_patches::Vector{SpecMor{BRT, BRET, RT, RET, MST1, MST2}}

  # fields used for caching
  inverse::SpecOpenMor{BRT, BRET, RT, RET, MST2, MST1}

  function SpecOpenMor(
      U::SpecOpen{BRT, BRET, RT, RET, MST1},
      V::SpecOpen{BRT, BRET, RT, RET, MST2},
      f::Vector{SpecMor{BRT, BRET, RT, RET, MST1, MST2}};
      check::Bool=true
    ) where {BRT, BRET, RT, RET, MST1, MST2}
    Y = parent(V)
    n = length(f)
    n == length(patches(U)) || error("number of patches does not coincide with the number of maps")
    if check
      for i in 1:n
	domain(f[i]) == patches(U)[i] || error("domain of definition of the map does not coincide with the patch")
	codomain(f[i]) == Y || error("codomain is not compatible")
      end
      for i in 1:n-1
	for j in i+1:n
	  A = intersect(domain(f[i]), domain(f[j]))
	  restrict(f[i], A, Y) == restrict(f[j], A, Y) || error("maps don't glue")
	end
      end
      I = ideal(localized_ring(OO(Y)), gens(V))
      for g in f
	one(localized_ring(OO(domain(g)))) in pullback(g)(I) + localized_modulus(OO(domain(g))) || error("image is not contained in the codomain")
      end
    end
    return new{BRT, BRET, RT, RET, MST1, MST2}(U, V, f)
  end
end

domain(f::SpecOpenMor) = f.domain
codomain(f::SpecOpenMor) = f.codomain
maps_on_patches(f::SpecOpenMor) = f.maps_on_patches
getindex(f::SpecOpenMor, i::Int) = maps_on_patches(f)[i]

@Markdown.doc """
maximal_extension(
  X::Spec{BRT, BRET, RT, RET, MST}, 
  Y::Spec{BRT, BRET, RT, RET, MST}, 
  f::AbstractAlgebra.Generic.Frac{RET}
) where {BRT, BRET, RT, RET, MST}

Given a rational map ϕ : X ---> Y ⊂ Spec 𝕜[y₁,…,yₙ] of affine schemes 
determined by ϕ*yⱼ = fⱼ = aⱼ//bⱼ, find the maximal open subset U⊂ X 
to which ϕ can be extended to a regular map g : U → Y and return g.
"""
function maximal_extension(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    Y::Spec{BRT, BRET, RT, RET, MST}, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {BRT, BRET, RT, RET, MST}
  U, g = maximal_extension(X, f)
  n = length(patches(U))
  maps = Vector{SpecMor{BRT, BRET, RT, RET, MST, MST}}()
  for i in 1:n
    push!(maps, SpecMor(patches(U)[i], Y, [representatives(a)[i] for a in g]))
  end
  return SpecOpenMor(U, SpecOpen(Y), maps)
end

function maximal_extension(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    Y::Spec{BRT, BRET, RT, RET, MST}, 
    f::Vector
  ) where {BRT, BRET, RT, RET, MST}
  return maximal_extension(X, Y, FractionField(base_ring(OO(X))).(f))
end

function restriction(
    f::SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2},
    U::SpecOpen{BRT, BRET, RT, RET, MST1},
    V::SpecOpen{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1, MST2}
  U == domain(f) || error("restriction not implemented")
  issubset(V, codomain(f))
  return SpecOpenMor(U, V, maps_on_patches(f))
end

