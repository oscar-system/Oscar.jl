export exceptional_prime_divisor

@doc raw"""
    underlying_morphism(bl::ToricBlowupMorphism)

Return the underlying toric morphism of a toric blowup. Access to other
attributes such as `domain`, `codomain`, `covering_morphism` are
executed via `underlying_morphism`.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 1, 1])
Toric blowup morphism

julia> Oscar.underlying_morphism(f)
Toric morphism
```
"""
underlying_morphism(bl::ToricBlowupMorphism) = bl.toric_morphism



@doc raw"""
    index_of_exceptional_ray(bl::ToricBlowupMorphism)

Return the index of the exceptional ray used in the construction of the toric
blowup.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 1, 1])
Toric blowup morphism

julia> index_of_exceptional_ray(f)
5
```
"""
index_of_exceptional_ray(bl::ToricBlowupMorphism) = bl.index_of_exceptional_ray


@doc raw"""
    minimal_supercone_coordinates_of_exceptional_ray(f::ToricBlowupMorphism) -> Vector{QQFieldElem}

Let $f\colon Y \to X$ be the toric blowup corresponding to a star
subdivision along a ray with minimal generator $v$.
This function returns the minimal supercone coordinate vector of $v$ in the fan of $X$.
See `?minimal_supercone_coordinates` for more details.

# Examples
```jldoctest
julia> X = affine_space(NormalToricVariety, 2)
Normal toric variety

julia> f = blow_up(X, [2, 3])
Toric blowup morphism

julia> minimal_supercone_coordinates_of_exceptional_ray(f)
2-element Vector{QQFieldElem}:
 2
 3
```
"""
@attr Vector{QQFieldElem} function minimal_supercone_coordinates_of_exceptional_ray(f::ToricBlowupMorphism)
  PF = polyhedral_fan(codomain(f))
  v = rays(domain(f))[index_of_exceptional_ray(f), :][1]
  v_ZZ = primitive_generator(v)
  return minimal_supercone_coordinates(PF, v_ZZ)
end

@doc raw"""
    exceptional_prime_divisor(bl::ToricBlowupMorphism)

Return the exceptional prime Weil divisor (as a toric divisor) of the
ray used to construct the toric blowup. Note that this divisor need not
be Cartier and this divisor need not coincide with the locus where the
morphism is not an isomorphism.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 2, 3])
Toric blowup morphism

julia> E = exceptional_prime_divisor(f)
Torus-invariant, prime divisor on a normal toric variety

julia> is_cartier(E)
false
```
"""
function exceptional_prime_divisor(bl::ToricBlowupMorphism)
  if !isdefined(bl, :exceptional_prime_divisor)
    X = domain(bl)
    S = cox_ring(X)
    x = gens(S)
    j = index_of_exceptional_ray(bl)
    help_list = [i == j ? 1 : 0 for i in 1:ngens(S)]
    td = toric_divisor(X, help_list)
    @assert is_prime(td) "exceptional prime divisor must be prime"
    bl.exceptional_prime_divisor = td
  end
  return bl.exceptional_prime_divisor
end

@doc raw"""
    center(bl::ToricBlowupMorphism) -> AbsIdealSheaf

Returns an ideal sheaf `I` such that the cosupport of `I` is the image
of the exceptional prime divisor.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 2, 3])
Toric blowup morphism

julia> center(f)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_3_1, x_2_1)
  2: Ideal (x_3_2, x_2_2)
  3: Ideal (1)
  4: Ideal (1)
```
"""
@attr AbsIdealSheaf function center(bl::ToricBlowupMorphism)
  X = domain(bl)
  S = cox_ring(X)
  # TODO: The current implementation is very slow.
  # Once ideal sheaves on nonsmooth normal toric varieties are
  # implemented, may replace the implementation with the following:
#   coords = minimal_supercone_coordinates_of_exceptional_ray(bl)
#   R = cox_ring(codomain(bl))
#   I = ideal(R, [gens(R)[i] for i in 1:ngens(R) if coords[i] > 0])
#   return ideal_sheaf(codomain(bl), I)
  x = gens(S)
  j = index_of_exceptional_ray(bl)
  I = ideal(S, x[j])
  II = ideal_sheaf(X, I)
  JJ = pushforward(bl, II)::IdealSheaf
  return JJ
end


#########################################################################
# Forwarding attributes of toric morphisms based on `underlying_morphism`
#########################################################################
lattice_homomorphism(bl::ToricBlowupMorphism) = lattice_homomorphism(underlying_morphism(bl))
morphism_on_torusinvariant_weil_divisor_group(bl::ToricBlowupMorphism) = morphism_on_torusinvariant_weil_divisor_group(underlying_morphism(bl))
morphism_on_torusinvariant_cartier_divisor_group(bl::ToricBlowupMorphism) = morphism_on_torusinvariant_cartier_divisor_group(underlying_morphism(bl))
morphism_on_class_group(bl::ToricBlowupMorphism) = morphism_on_class_group(underlying_morphism(bl))
morphism_on_picard_group(bl::ToricBlowupMorphism) = morphism_on_picard_group(underlying_morphism(bl))



#=
For the future, if traits seem better to forward methods as the ones immediately above (lattice_homomorphism, morphism_on_torusinvariant_weil_divisor_group, etc.):

########################################################################
# Enabling the HasToricSubObjectTrait for ToricBlowupMorphism        #
########################################################################
HasToricSubObjectTrait(::Type{T}) where {T<:ToricBlowupMorphism} = HasToricSubObjectTrait{ToricMorphism}()

toric_sub_object(bl::ToricBlowupMorphism) = underlying_morphism(bl)

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
