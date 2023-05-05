@doc raw"""
    global_tate_model_over_projective_space()

This method constructs a global Tate model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> global_tate_model_over_projective_space()
Global Tate model over a concrete base
```
"""
global_tate_model_over_projective_space() = global_tate_model(projective_space(NormalToricVariety,3))


@doc raw"""
    global_weierstrass_model_over_projective_space()

This method constructs a global Weierstrass model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> global_weierstrass_model_over_projective_space()
Global Weierstrass model over a concrete base
```
"""
global_weierstrass_model_over_projective_space() = global_weierstrass_model(projective_space(NormalToricVariety,3))


@doc raw"""
    su5_tate_model_over_arbitrary_3d_base()

Return the SU(5) Tate model over an arbitrary
3-dimensional base space. For more details see
e.g. [Wei18](@cite) and references therein.

```jldoctest
julia> tm = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> v = underlying_toric_variety(ambient_space(tm))
Normal toric variety

julia> a10,a21,a32,a43,a65,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = blow_up(v,I)
Normal toric variety

julia> cox_ring(v2)
Multivariate Polynomial Ring in 10 variables a10, a21, a32, a43, ..., e over Rational Field graded by 
  a10 -> [0 0]
  a21 -> [0 0]
  a32 -> [0 0]
  a43 -> [0 0]
  a65 -> [0 0]
  w -> [1 0]
  x -> [1 2]
  y -> [1 3]
  z -> [0 1]
  e -> [-1 0]
```
"""
function su5_tate_model_over_arbitrary_3d_base()
    auxiliary_base_ring, (a10, a21, a32, a43, a65, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];
    a1 = a10;
    a2 = a21 * w;
    a3 = a32 * w^2;
    a4 = a43 * w^3;
    a6 = a65 * w^5;
    ais = [a1, a2, a3, a4, a6];
    return global_tate_model(ais, auxiliary_base_ring, 3)
end

@doc raw"""
    su5_weierstrass_model_over_arbitrary_3d_base()

Return the SU(5) Weierstrass model over an arbitrary
3-dimensional base space. For more details see
e.g. [Wei18](@cite) and references therein.

```jldoctest
julia> tm = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base
```
"""
su5_weierstrass_model_over_arbitrary_3d_base() = global_weierstrass_model(su5_tate_model_over_arbitrary_3d_base())
