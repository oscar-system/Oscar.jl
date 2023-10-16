@doc raw"""
    total_transform(f::AbsSimpleBlowdownMorphism, II::IdealSheaf)

Computes the total transform of an ideal sheaf along a blowdown morphism.

In particular, this applies in the toric setting. However, note that
currently (October 2023), ideal sheaves are only supported on smooth
toric varieties.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> bl = blow_up(P2, [1, 1])
Toric blowdown morphism

julia> S = cox_ring(P2);

julia> x, y, z = gens(S);

julia> I = ideal_sheaf(P2, ideal([x*y]))
Sheaf of ideals
  on normal, smooth toric variety
with restrictions
  1: ideal(x_1_1*x_2_1)
  2: ideal(x_2_2)
  3: ideal(x_1_3)

julia> total_transform(bl, I)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: ideal(x_1_1*x_2_1^2)
  2: ideal(x_1_2^2*x_2_2)
  3: ideal(x_2_3)
  4: ideal(x_1_4)
```
"""
function total_transform(f::AbsSimpleBlowdownMorphism, II::IdealSheaf)
  return pullback(f, II)
end

function total_transform(f::AbsBlowdownMorphism, II::IdealSheaf)
  return pullback(f, II)
end
