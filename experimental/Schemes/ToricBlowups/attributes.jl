@doc raw"""
    underlying_morphism(bl::ToricBlowdownMorphism)

Return the underlying toric morphism of a toric blowdown morphism.
Access to other attributes such as `domain`, `codomain`, `covering_morphism`
are executed via `underlying_morphism`.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> blow_down_morphism = blow_up(P3, [0, 1, 1])
Toric blowdown morphism

julia> Oscar.underlying_morphism(blow_down_morphism)
Toric morphism
```
"""
underlying_morphism(bl::ToricBlowdownMorphism) = bl.toric_morphism



@doc raw"""
    index_of_new_ray(bl::ToricBlowdownMorphism)

Return the index of the new ray used in the construction of the toric blowdown morphism.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> blow_down_morphism = blow_up(P3, [0, 1, 1])
Toric blowdown morphism

julia> index_of_new_ray(blow_down_morphism)
5
```
"""
index_of_new_ray(bl::ToricBlowdownMorphism) = bl.index_of_new_ray



@doc raw"""
    center(bl::ToricBlowdownMorphism)

Return the center of the toric blowdown morphism as ideal sheaf.

Currently (October 2023), ideal sheaves are only supported for
smooth toric varieties. Hence, there will be instances in which
the construction of a blowdown morphism succeeds, but the center
of the blowup, as represented by an ideal sheaf, cannot yet be
computed.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> blow_down_morphism = blow_up(P3, [0, 1, 1])
Toric blowdown morphism

julia> center(blow_down_morphism)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_3_1, x_2_1)
  2: Ideal (x_3_2, x_2_2)
  3: Ideal (1)
  4: Ideal (1)
```
"""
function center(bl::ToricBlowdownMorphism)
  if !isdefined(bl, :center)
    # TODO: The implementation below is highly inefficient. Improve it if you know how.
    X = domain(bl)
    S = cox_ring(X)
    x = gens(S)
    j = index_of_new_ray(bl)
    I = ideal(S, x[j])
    II = ideal_sheaf(X, I)
    JJ = pushforward(bl, II)::IdealSheaf
    bl.center = JJ
  end
  return bl.center
end



@doc raw"""
    exceptional_divisor(bl::ToricBlowdownMorphism)

Return the exceptional divisor (as toric divisor)  of the toric blowdown morphism.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> blow_down_morphism = blow_up(P3, [0, 1, 1])
Toric blowdown morphism

julia> exceptional_divisor(blow_down_morphism)
Torus-invariant, cartier, prime divisor on a normal toric variety
```
"""
function exceptional_divisor(bl::ToricBlowdownMorphism)
  if !isdefined(bl, :exceptional_divisor)
    X = domain(bl)
    S = cox_ring(X)
    x = gens(S)
    j = index_of_new_ray(bl)
    help_list = [i == j ? 1 : 0 for i in 1:ngens(S)]
    td = toric_divisor(X, help_list)
    @assert is_cartier(td) "exceptional divisor must be Cartier"
    @assert is_prime(td) "exceptional divisor must be prime"
    bl.exceptional_divisor = td
  end
  return bl.exceptional_divisor
end



#########################################################################
# Forwarding attributes of toric morphisms based on `underlying_morphism`
#########################################################################
grid_morphism(bl::ToricBlowdownMorphism) = grid_morphism(underlying_morphism(bl))
morphism_on_torusinvariant_weil_divisor_group(bl::ToricBlowdownMorphism) = morphism_on_torusinvariant_weil_divisor_group(underlying_morphism(bl))
morphism_on_torusinvariant_cartier_divisor_group(bl::ToricBlowdownMorphism) = morphism_on_torusinvariant_cartier_divisor_group(underlying_morphism(bl))
morphism_on_class_group(bl::ToricBlowdownMorphism) = morphism_on_class_group(underlying_morphism(bl))
morphism_on_picard_group(bl::ToricBlowdownMorphism) = morphism_on_picard_group(underlying_morphism(bl))



#=
For the future, if traits seem better to forward methods as the ones immediately above (grid_morphism, morphism_on_torusinvariant_weil_divisor_group, etc.):

########################################################################
# Enabling the HasToricSubObjectTrait for ToricBlowdownMorphism        #
########################################################################
HasToricSubObjectTrait(::Type{T}) where {T<:ToricBlowdownMorphism} = HasToricSubObjectTrait{ToricMorphism}()

toric_sub_object(bl::ToricBlowdownMorphism) = underlying_morphism(bl)

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
