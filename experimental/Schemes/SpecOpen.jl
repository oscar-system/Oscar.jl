export SpecOpen, ambient, gens, complement, npatches, patches, intersections, name, intersect, issubset, closure, find_non_zero_divisor, is_non_zero_divisor, is_dense

export StructureSheafRing, variety, domain, OO

export StructureSheafElem, domain, representatives, patches, restrict, npatches

export SpecOpenMor, maps_on_patches, restriction, identity_map, preimage, generic_fractions, pullback, maximal_extension

@Markdown.doc """
    SpecOpen{BRT, BRET, RT, RET, MST} <: Scheme{BRT, BRET}

Zariski open subset ``U`` of an affine scheme ``X = Spec(R)``. 
This stores a list of generators ``f₁,…,fᵣ`` of an ideal 
``I`` defining the complement ``Z = X ∖ U``. 
The variety ``X`` is referred to as the *ambient variety* and 
the list ``f₁,…,fᵣ`` as the *generators* for ``U``.
"""
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

@Markdown.doc """
    ambient(U::SpecOpen)

Returns the ambient scheme ``X`` of a Zariski open subset ``U ⊂ X``.
"""
ambient(U::SpecOpen) = U.X

@Markdown.doc """
    npatches(U::SpecOpen)

Returns the number of generators stored for describing the complement of ``U``.
"""
npatches(U::SpecOpen) = length(U.gens)

@Markdown.doc """
    gens(U::SpecOpen) = U.gens

Returns the generators ``[f₁,…,fᵣ]`` stored for the description 
of the complement of ``U``.
"""
gens(U::SpecOpen) = U.gens

@Markdown.doc """
    getindex(U::SpecOpen, i::Int) = patches(U)[i]

Returns the hypersurface complement of ``fᵢ`` in the 
ambient scheme ``X`` of ``U`` where ``f₁,…,fᵣ`` are 
the generators stored for the description of the complement 
of ``U``.
"""
getindex(U::SpecOpen, i::Int) = patches(U)[i]

function Base.show(io::IO, U::SpecOpen)
  if isdefined(U, :name) 
    print(io, name(U))
    return
  end
  print(io, "complement of zero locus of $(gens(U)) in $(ambient(U))")
end

@Markdown.doc """
    function SpecOpen(
      X::Spec{BRT, BRET, RT, RET, MST},
      I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
    ) where {BRT, BRET, RT, RET, MST}

Returns the complement of the zero locus of ``I`` in ``X``.
"""
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
    I = radical(saturated_ideal(ideal(localized_ring(OO(ambient(U))), gens(U))))
    U.complement = subscheme(ambient(U), I)
  end
  return U.complement
end

@Markdown.doc """
    patches(U::SpecOpen)

Returns a list of principal affine open subschemes covering ``U``.
"""
function patches(U::SpecOpen)
  if !isdefined(U, :patches)
    X = ambient(U)
    V = Vector{typeof(X)}()
    for f in gens(U)
      push!(V, hypersurface_complement(X, f))
    end
    U.patches = V
  end
  return U.patches
end

@Markdown.doc """
    intersections(U::SpecOpen)

Returns a list of pairwise intersections of the 
principal open subschemes covering ``U``.
"""
function intersections(U::SpecOpen)
  if !isdefined(U, :intersections)
    X = ambient(U)
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
  return "open subset of $(ambient(U))"
end

function intersect(
    Y::Spec{BRT, BRET, RT, RET, MST}, 
    U::SpecOpen{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  X = ambient(U)
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

function intersect(
    U::SpecOpen{BRT, BRET, RT, RET, MST},
    V::SpecOpen{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  X = ambient(U) 
  X == ambient(V) || error("method not implemented")
  return SpecOpen(X, [a*b for a in gens(U) for b in gens(V)])
end

function Base.union(U::T, V::T) where {T<:SpecOpen}
  ambient(U) == ambient(V) || error("the two open sets do not lay in the same ambient variety")
  return SpecOpen(ambient(U), vcat(gens(U), gens(V)))
end

function issubset(
    Y::Spec{BRT, BRET, RT, RET, MST}, 
    U::SpecOpen{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  return one(localized_ring(OO(Y))) in ideal(OO(Y), gens(U))
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
  base_ring(OO(ambient(U))) == base_ring(OO(ambient(V))) || return false
  return issubset(complement(intersect(V, ambient(U))), complement(U))
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

@Markdown.doc """
    closure(U::SpecOpen)

Computes the Zariski closure of an open set ``U ⊂ X`` 
with ``X`` is the affine ambient scheme of ``U``.
"""
function closure(U::SpecOpen)
  X = ambient(U)
  R = base_ring(OO(X))
  I = saturated_ideal(localized_modulus(OO(X)))
  I = saturation(I, ideal(R, gens(U)))
  return subscheme(X, I)
end

@Markdown.doc """
    closure(
      U::SpecOpen{BRT, BRET, RT, RET, MST},
      Y::Spec{BRT, BRET, RT, RET, MST} 
    ) where {BRT, BRET, RT, RET, MST}

Computes the closure of ``U ⊂ Y``.
"""
function closure(
    U::SpecOpen{BRT, BRET, RT, RET, MST},
    Y::Spec{BRT, BRET, RT, RET, MST} 
  ) where {BRT, BRET, RT, RET, MST}
  issubset(U, Y) || error("the first set is not contained in the second")
  X = closure(U)
  return intersect(X, Y)
end


@Markdown.doc """
StructureSheafRing{BRT, BRET, RT, RET, MST}

The ring of regular functions ``𝒪(X, U)`` on an open subset ``U`` of an 
affine scheme ``X``.
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

@Markdown.doc """
    variety(R::StructureSheafRing)

The ring ``R = 𝒪(X, U)`` belongs to a sheaf of rings ``𝒪(X, -)`` and this returns 
the variety ``X`` on which ``𝒪`` is defined.
"""
variety(R::StructureSheafRing) = R.variety

@Markdown.doc """
    domain(R::StructureSheafRing)

For a ring ``R = 𝒪(X, U)`` this returns ``U``.
"""
domain(R::StructureSheafRing) = R.domain

OO(U::SpecOpen) = StructureSheafRing(ambient(U), U)
OO(X::Spec{BRT, BRET, RT, RET, MST}, U::SpecOpen{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = StructureSheafRing(X, U)

function ==(R::T, S::T) where {T<:StructureSheafRing} 
  variety(R) == variety(S) || return false
  domain(S) == domain(R) || return false
  return true
end

elem_type(R::StructureSheafRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}= Type{StructureSheafElem{BRT, BRET, RT, RET, MST}}

@Markdown.doc """
StructureSheafElem{BRT, BRET, RT, RET, MST}

An element ``f ∈ 𝒪(X, U)`` of the ring of regular functions on 
an open set ``U`` of an affine variety ``X``.
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

variety(f::StructureSheafElem) = ambient(domain(f))
domain(f::StructureSheafElem) = f.domain
representatives(f::StructureSheafElem) = f.representatives
patches(f::StructureSheafElem) = patches(domain(f))
npatches(f::StructureSheafElem) = length(f.representatives)
getindex(f::StructureSheafElem, i::Int) = getindex(representatives(f), i)
ambient(f::StructureSheafElem) = OO(variety(f), domain(f))

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
  VU = [intersect(V, U) for U in patches(domain(f))]
  g = [OO(VU[i])(f[i]) for i in 1:length(VU)]
  l = write_as_linear_combination(one(OO(V)), OO(V).(lifted_denominator.(g)))
  return dot(l, OO(V).(lifted_numerator.(g)))
end

@Markdown.doc """
    maximal_extension(
      X::Spec{BRT, BRET, RT, RET, MST}, 
      f::AbstractAlgebra.Generic.Frac{RET}
    ) where {BRT, BRET, RT, RET, MST}

Returns the maximal extension of the restriction of ``f`` 
to a rational function on ``X`` on a maximal domain of 
definition ``U ⊂ X``. 
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

Returns the extension of the restriction of the ``fᵢ`` as a 
set of rational functions on ``X`` as *regular* functions to a 
common maximal domain of definition ``U ⊂ X``.
"""
function maximal_extension(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {BRT, BRET, RT, RET, MST}
  if length(f) == 0
    return SpecOpen(X), Vector{StructureSheafElem{BRT, BRET, RT, RET, MST}}()
  end
  R = base_ring(parent(f[1]))
  for a in f
    R == base_ring(parent(a)) || error("fractions do not belong to the same ring")
  end
  R == base_ring(OO(X)) || error("fractions do not belong to the base ring of the scheme")
  a = numerator.(f)
  b = denominator.(f)
  W = localized_ring(OO(X))
  I = ideal(W, one(W))
  for p in f
    I = intersect(quotient(ideal(W, denominator(p)) + localized_modulus(OO(X)), ideal(W, numerator(p))), I)
  end
  U = SpecOpen(X, I)
  return U, [StructureSheafElem(U, [OO(V)(a) for V in patches(U)]) for a in f]
end

@Markdown.doc """
    SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2}

Morphisms ``f : U → V`` of open sets ``U ⊂ X`` and ``V ⊂ Y`` of affine schemes.
These are stored as morphisms ``fᵢ: Uᵢ→ Y`` on the affine patches 
``Uᵢ`` of ``U``.
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
    Y = ambient(V)
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

function SpecOpenMor(U::SpecOpenType, V::SpecOpenType, f::Vector) where {SpecOpenType<:SpecOpen}
  Y = ambient(V)
  return SpecOpenMor(U, V, [SpecMor(W, Y, f) for W in patches(U)]) 
end

function inclusion_morphism(U::SpecOpenType, V::SpecOpenType) where {SpecOpenType<:SpecOpen}
  X = ambient(U)
  ambient(V) == X || error("method not implemented")
  return SpecOpenMor(U, V, gens(base_ring(OO(X))))
end


@Markdown.doc """
    compose(f::T, g::T) where {T<:SpecOpenMor}

computes the composition of two morphisms 

      f    g
    U →  V →  W
    ∩    ∩    ∩
    X    Y    Z

of open sets of affine varieties.
"""
function compose(f::T, g::T) where {T<:SpecOpenMor}
  U = domain(f)
  Cf = codomain(f)
  V = domain(g)
  issubset(Cf, V) || error("maps are not compatible")
  W = codomain(g)
  X = ambient(U)
  Y = ambient(V)
  Z = ambient(W)
  g_maps = maps_on_patches(g)
  f_maps = maps_on_patches(f)
  #####################################################################
  # The method proceeds as follows.
  # We extract the affine open sets 
  #   Uᵢⱼ = Uᵢ∩ f⁻¹(Vⱼ)  
  # where Uᵢ are the patches of U and Vⱼ those of V.
  # Then we can pass to the restrictions 
  #   fᵢⱼ : Uᵢⱼ → Vⱼ 
  # and the compositions 
  #   gᵢⱼ :=  gⱼ ∘ fᵢⱼ: Uᵢⱼ → Z.
  # For any variable z for Z we can write the pullbacks gᵢⱼ*z ∈ 𝒪(Uᵢⱼ). 
  # These functions necessarily coincide on the overlaps of the Uᵢⱼ in 
  # Uᵢ, so we can extend them to a global function on Uᵢ. 
  # From these global functions, we can then assemble a morphism 
  # on each one of the affine schemes Uᵢ which glue together to a 
  # SpecOpenMor over these patches.
  #####################################################################
  m = length(gens(U))
  n = length(gens(V))
  result_maps = Vector{typeof(f_maps[1])}()
  for i in 1:m
    U_i = patches(U)[i]
    f_i = f_maps[i]
    z_i = Vector{elem_type(OO(U_i))}() # the pullbacks of the coordinate functions of Z on Uᵢ
    z_ij = Vector{Vector{elem_type(OO(U_i))}}()
    for j in 1:n
      V_j = patches(V)[j]
      g_j = g_maps[j]
      U_ij = intersect(U_i, preimage(f_i, V_j))
      g_ij = compose(restrict(f_i, U_ij, V_j), g_j)
      push!(z_ij, [pullback(g_ij)(z) for z in gens(OO(Z))])
      ### This is using the well known trick that if aᵢ// dᵢ coincide pairwise
      # wherever dᵢ and dⱼ are defined, and 1 = ∑ᵢ λᵢ⋅ dᵢ, then 
      # ∑ᵢλᵢ⋅ aᵢ = aᵢ// dᵢ on all open sets {dᵢ≠ 0}. 
    end
    for k in 1:length(gens(OO(Z)))
      l = write_as_linear_combination(one(OO(U_i)), lifted_denominator.([z_ij[j][k] for j in 1:n]))
      push!(z_i, dot(l, OO(U_i).(lifted_numerator.([z_ij[j][k] for j in 1:n]))))
    end
    push!(result_maps, SpecMor(U_i, Z, MPolyQuoLocalizedRingHom(OO(Z), OO(U_i), z_i)))
    #push!(result_maps, SpecMor(U_i, Z, z_i))
  end
  return SpecOpenMor(U, W, result_maps)
end

#function restriction(
#   f::SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2}, 
#   U::SpecOpen{BRT, BRET, RT, RET, MST1},
#   V::SpecOpen{BRT, BRET, RT, RET, MST2}
# ) where {BRT, BRET, RT, RET, MST1, MST2}
# g = compose(inclusion_morphism(U, domain(f)), f)
# return SpecOpenMor(U, V, maps_on_patches(g))
#end

function pullback(f::SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2}, a::RET) where {BRT, BRET, RT, RET, MST1, MST2}
  U = domain(f)
  X = ambient(U)
  V = codomain(f)
  Y = ambient(V)
  R = base_ring(OO(Y))
  parent(a) == R || error("element does not belong to the correct ring")
  pb_a = [pullback(f[i])(a) for i in 1:npatches(U)]
  return StructureSheafElem(U, pb_a)
end

@Markdown.doc """
    maximal_extension(
      X::Spec{BRT, BRET, RT, RET, MST}, 
      Y::Spec{BRT, BRET, RT, RET, MST}, 
      f::AbstractAlgebra.Generic.Frac{RET}
    ) where {BRT, BRET, RT, RET, MST}

Given a rational map ``ϕ : X ---> Y ⊂ Spec 𝕜[y₁,…,yₙ]`` of affine schemes 
determined by ``ϕ*(yⱼ) = fⱼ = aⱼ/bⱼ``, find the maximal open subset ``U⊂ X`` 
to which ``ϕ`` can be extended to a regular map ``g : U → Y`` and return ``g``.
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
  inc = SpecOpenMor(U, domain(f), [SpecMor(W, ambient(domain(f)), gens(localized_ring(OO(W)))) for W in patches(U)])
  help_map = compose(inc, f)
  return SpecOpenMor(U, V, maps_on_patches(help_map))
end

identity_map(U::SpecOpen) = SpecOpenMor(U, U, [SpecMor(V, ambient(U), gens(localized_ring(OO(V)))) for V in patches(U)])

function ==(f::T, g::T) where {T<:SpecOpenMor} 
  domain(f) == domain(g) || return false
  codomain(f) == codomain(g) || return false
  Y = ambient(codomain(f))
  m = length(patches(domain(f)))
  n = length(patches(domain(g)))
  for i in 1:m
    for j in 1:n
      restrict(f[i], intersect(domain(f)[i], domain(g)[i]), Y) == 
      restrict(g[i], intersect(domain(f)[i], domain(g)[i]), Y) || return false
    end
  end
  return true
end

function preimage(f::SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2},
    Z::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1, MST2}
  U = domain(f) 
  X = ambient(U)
  n = length(patches(U))
  W = localized_ring(OO(X))
  I = ideal(W, one(W))
  for i in 1:n 
    I = intersect(I, localized_modulus(OO(closure(preimage(f[i], Z), X))))
  end
  fZbar = subscheme(X, I)
  return SpecOpen(fZbar, gens(U))
end

function preimage(f::SpecOpenMor{BRT, BRET, RT, RET, MST1, MST2},
    V::SpecOpen{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1, MST2}
  U = domain(f)
  X = ambient(U)
  R = base_ring(OO(X))
  I = ideal(R, one(R))
  for i in 1:npatches(U)
    I = intersect(I, saturated_ideal(ideal(OO(U[i]), pullback(f[i]).(gens(V)))))
  end
  return intersect(U, SpecOpen(X, gens(I)))
end

function is_non_zero_divisor(f::RET, U::SpecOpen{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET}
  for V in patches(U)
    zero_ideal = ideal(OO(V), [zero(OO(V))])
    zero_ideal == quotient(zero_ideal, ideal(OO(V), [f])) || return false
  end
  return true
end

function find_non_zero_divisor(U::SpecOpen)
  n = length(gens(U))
  X = ambient(U)
  kk = coefficient_ring(base_ring(OO(X)))
  d = dot([rand(kk, 0:100) for i in 1:n], gens(U))
  while !is_non_zero_divisor(d, U)
    d = dot([rand(kk, 0:100) for i in 1:n], gens(U))
  end
  return d
end

function generic_fractions(f::SpecOpenMor)
  U = domain(f)
  X = ambient(U)
  V = codomain(f)
  Y = ambient(V)
  d = find_non_zero_divisor(U)
  V = hypersurface_complement(X, d)
  result = fraction.([restrict(pullback(f, y), V) for y in gens(base_ring(OO(Y)))])
end

function is_dense(U::SpecOpen)
  X = ambient(U)
  I = [localized_modulus(OO(closure(V, X))) for V in patches(U)]
  J = ideal(OO(X), [one(OO(X))])
  for i in I
    J = intersect(J, i)
  end
  return J == localized_modulus(OO(X))
end
