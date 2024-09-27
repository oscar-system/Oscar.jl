export center_data
export center_unnormalized
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
    index_of_new_ray(bl::ToricBlowupMorphism)

Return the index of the new ray used in the construction of the toric
blowup.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 1, 1])
Toric blowup morphism

julia> index_of_new_ray(f)
5
```
"""
index_of_new_ray(bl::ToricBlowupMorphism) = bl.index_of_new_ray


@doc raw"""
    center_data(bl::ToricBlowupMorphism)

Returns the ideal, ideal sheaf or ray that was used to construct the
morphism.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 2, 3])
Toric blowup morphism

julia> center_data(f)
3-element Vector{Int64}:
 0
 2
 3
```
"""
function center_data(bl::ToricBlowupMorphism)
  return bl.center_data
end


@doc raw"""
    center_unnormalized(bl::ToricBlowupMorphism)

Returns an ideal sheaf `I` such that the normalization of the blowup
along `I` gives the morphism `bl`.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> f = blow_up(P3, [0, 2, 3])
Toric blowup morphism

julia> center_unnormalized(f)
Sheaf of ideals
  on normal, smooth toric variety
with restrictions
  1: Ideal (x_2_1^2, x_3_1^3)
  2: Ideal (x_2_2^2, x_3_2^3)
  3: Ideal (1)
  4: Ideal (1)
```
"""
function center_unnormalized(bl::ToricBlowupMorphism)
  if !isdefined(bl, :center_unnormalized)
    # TODO: The implementation below is highly inefficient. Improve it if you know how.
    X = domain(bl)
    S = cox_ring(X)
    x = gens(S)
    j = index_of_new_ray(bl)
    I = ideal(S, x[j])
    II = ideal_sheaf(X, I)
    JJ = pushforward(bl, II)::IdealSheaf
    bl.center_unnormalized = JJ
  end
  return bl.center_unnormalized
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
    j = index_of_new_ray(bl)
    help_list = [i == j ? 1 : 0 for i in 1:ngens(S)]
    td = toric_divisor(X, help_list)
    @assert is_prime(td) "exceptional prime divisor must be prime"
    bl.exceptional_prime_divisor = td
  end
  return bl.exceptional_prime_divisor
end



#########################################################################
# Forwarding attributes of toric morphisms based on `underlying_morphism`
#########################################################################
grid_morphism(bl::ToricBlowupMorphism) = grid_morphism(underlying_morphism(bl))
morphism_on_torusinvariant_weil_divisor_group(bl::ToricBlowupMorphism) = morphism_on_torusinvariant_weil_divisor_group(underlying_morphism(bl))
morphism_on_torusinvariant_cartier_divisor_group(bl::ToricBlowupMorphism) = morphism_on_torusinvariant_cartier_divisor_group(underlying_morphism(bl))
morphism_on_class_group(bl::ToricBlowupMorphism) = morphism_on_class_group(underlying_morphism(bl))
morphism_on_picard_group(bl::ToricBlowupMorphism) = morphism_on_picard_group(underlying_morphism(bl))



#=
For the future, if traits seem better to forward methods as the ones immediately above (grid_morphism, morphism_on_torusinvariant_weil_divisor_group, etc.):

########################################################################
# Enabling the HasToricSubObjectTrait for ToricBlowupMorphism        #
########################################################################
HasToricSubObjectTrait(::Type{T}) where {T<:ToricBlowupMorphism} = HasToricSubObjectTrait{ToricMorphism}()

toric_sub_object(bl::ToricBlowupMorphism) = underlying_morphism(bl)

########################################################################
# Enabling the functionality in general                                #
########################################################################
function grid_morphism(phi::Any)
  return _grid_morphism(HasToricSubObjectTrait(phi), phi)
end

function _grid_morphism(::HasToricSubObject, phi)
  return grid_morphism(toric_sub_object(phi))
end
=#
