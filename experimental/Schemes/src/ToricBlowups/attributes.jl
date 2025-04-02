export exceptional_prime_divisor

# Return the underlying toric morphism of a toric blowup. Access to other
# attributes such as `domain`, `codomain`, `covering_morphism` are
# executed via `underlying_morphism`.
# Example:
# ```jldoctest
# julia> X = projective_space(NormalToricVariety, 3)
# Normal toric variety
#
# julia> phi = blow_up(X, [0, 1, 1])
# Toric blowup morphism
#
# julia> Oscar.underlying_morphism(phi)
# Toric morphism
# ```
# """
underlying_morphism(phi::ToricBlowupMorphism) = phi.toric_morphism

@doc raw"""
    index_of_exceptional_ray(phi::ToricBlowupMorphism) -> Int

Return the index of the exceptional ray used in the construction of the toric
blowup.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> phi = blow_up(X, [0, 1, 1])
Toric blowup morphism

julia> index_of_exceptional_ray(phi)
5
```
"""
index_of_exceptional_ray(phi::ToricBlowupMorphism) = phi.index_of_exceptional_ray

@doc raw"""
    minimal_supercone_coordinates_of_exceptional_ray(phi::ToricBlowupMorphism) -> Vector{QQFieldElem}

Let $\varphi\colon Y \to X$ be the toric blowup corresponding to a star
subdivision along a primitive vector $r$ in the support of the fan of $X$.
This function returns the minimal supercone coordinate vector of $r$
(the output of `minimal_supercone_coordinates(polyhedral_fan(X), r)`).

# Examples
```jldoctest
julia> X = affine_space(NormalToricVariety, 2)
Normal toric variety

julia> phi = blow_up(X, [2, 3])
Toric blowup morphism

julia> minimal_supercone_coordinates_of_exceptional_ray(phi)
2-element Vector{QQFieldElem}:
 2
 3
```
"""
@attr Vector{QQFieldElem} function minimal_supercone_coordinates_of_exceptional_ray(phi::ToricBlowupMorphism)
  fan = polyhedral_fan(codomain(phi))
  r = rays(domain(phi))[index_of_exceptional_ray(phi), :][1]
  r_ZZ = primitive_generator(r)
  return minimal_supercone_coordinates(fan, r_ZZ)
end

@doc raw"""
    exceptional_prime_divisor(phi::ToricBlowupMorphism) -> ToricDivisor

Return the exceptional prime Weil divisor (as a toric divisor) of the
ray used to construct the toric blowup. Note that this divisor need not
be Cartier and this divisor need not coincide with the locus where the
morphism is not an isomorphism.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> phi = blow_up(X, [0, 2, 3])
Toric blowup morphism

julia> E = exceptional_prime_divisor(phi)
Torus-invariant, prime divisor on a normal toric variety

julia> is_cartier(E)
false
```
"""
function exceptional_prime_divisor(phi::ToricBlowupMorphism)
  if !isdefined(phi, :exceptional_prime_divisor)
    X = domain(phi)
    S = cox_ring(X)
    x = gens(S)
    j = index_of_exceptional_ray(phi)
    help_list = [i == j ? 1 : 0 for i in 1:ngens(S)]
    E = toric_divisor(X, help_list)
    @req is_prime(E) "exceptional prime divisor must be prime"
    phi.exceptional_prime_divisor = E
  end
  return phi.exceptional_prime_divisor
end

@doc raw"""
    center(phi::ToricBlowupMorphism) -> AbsIdealSheaf

Return an ideal sheaf $\mathcal{I}$ such that the cosupport of
$\mathcal{I}$ is the image of the exceptional prime divisor.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> phi = blow_up(X, [0, 2, 3])
Toric blowup morphism

julia> center(phi)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_3_1, x_2_1)
  2: Ideal (x_3_2, x_2_2)
  3: Ideal (1)
  4: Ideal (1)
```
"""
@attr AbsIdealSheaf function center(phi::ToricBlowupMorphism)
  X = domain(phi)
  S = cox_ring(X)
  # TODO: The current implementation is very slow.
  # Once ideal sheaves on nonsmooth normal toric varieties are
  # implemented, may replace the implementation with the following:
#   coords = minimal_supercone_coordinates_of_exceptional_ray(bl)
#   R = cox_ring(codomain(bl))
#   I = ideal(R, [gens(R)[i] for i in 1:ngens(R) if coords[i] > 0])
#   return ideal_sheaf(codomain(bl), I)
  x = gens(S)
  j = index_of_exceptional_ray(phi)
  I = ideal(S, x[j])
  II = ideal_sheaf(X, I)
  JJ = pushforward(phi, II)::IdealSheaf
  return JJ
end


#########################################################################
# Forwarding attributes of toric morphisms based on `underlying_morphism`
#########################################################################
lattice_homomorphism(phi::ToricBlowupMorphism) = lattice_homomorphism(underlying_morphism(phi))
morphism_on_torusinvariant_weil_divisor_group(phi::ToricBlowupMorphism) = morphism_on_torusinvariant_weil_divisor_group(underlying_morphism(phi))
morphism_on_torusinvariant_cartier_divisor_group(phi::ToricBlowupMorphism) = morphism_on_torusinvariant_cartier_divisor_group(underlying_morphism(phi))
morphism_on_class_group(phi::ToricBlowupMorphism) = morphism_on_class_group(underlying_morphism(phi))
morphism_on_picard_group(phi::ToricBlowupMorphism) = morphism_on_picard_group(underlying_morphism(phi))



#=
For the future, if traits seem better to forward methods as the ones immediately above (lattice_homomorphism, morphism_on_torusinvariant_weil_divisor_group, etc.):

########################################################################
# Enabling the HasToricSubObjectTrait for ToricBlowupMorphism        #
########################################################################
HasToricSubObjectTrait(::Type{T}) where {T<:ToricBlowupMorphism} = HasToricSubObjectTrait{ToricMorphism}()

toric_sub_object(phi::ToricBlowupMorphism) = underlying_morphism(phi)

########################################################################
# Enabling the functionality in general                                #
########################################################################
function lattice_homomorphism(phi::Any)
  return _lattice_homomorphism(HasToricSubObjectTrait(phi), phi)
end

function _lattice_homomorphism(::HasToricSubObject, phi)
  return lattice_homomorphism(toric_sub_object(phi))
end
=#
