#####################################################
# 1. Models over projective space
#####################################################

@doc raw"""
    global_tate_model_over_projective_space(d::Int)

This method constructs a global Tate model over the projective space.

# Examples
```jldoctest
julia> global_tate_model_over_projective_space(3)
Global Tate model over a concrete base
```
"""
global_tate_model_over_projective_space(d::Int) = global_tate_model(projective_space(NormalToricVariety, d); completeness_check = false)


@doc raw"""
    weierstrass_model_over_projective_space(d::Int)

This method constructs a Weierstrass model over the projective space.

# Examples
```jldoctest
julia> weierstrass_model_over_projective_space(3)
Weierstrass model over a concrete base
```
"""
weierstrass_model_over_projective_space(d::Int) = weierstrass_model(projective_space(NormalToricVariety, d); completeness_check = false)


@doc raw"""
    hypersurface_model_over_projective_space(d::Int)

This method constructs a hypersurface model over the projective space.

# Examples
```jldoctest
julia> hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base
```
"""
hypersurface_model_over_projective_space(d::Int) = hypersurface_model(projective_space(NormalToricVariety, d); completeness_check = false)


#####################################################
# 2. Models over hirzebruch surfaces
#####################################################

@doc raw"""
    global_tate_model_over_hirzebruch_surface(r::Int)

This method constructs a global Tate model over a Hirzebruch surface.

# Examples
```jldoctest
julia> global_tate_model_over_hirzebruch_surface(1)
Global Tate model over a concrete base
```
"""
global_tate_model_over_hirzebruch_surface(r::Int) = global_tate_model(hirzebruch_surface(NormalToricVariety, r); completeness_check = false)


@doc raw"""
    weierstrass_model_over_hirzebruch_surface(r::Int)

This method constructs a Weierstrass model over a Hirzebruch surface.

# Examples
```jldoctest
julia> weierstrass_model_over_hirzebruch_surface(1)
Weierstrass model over a concrete base
```
"""
weierstrass_model_over_hirzebruch_surface(r::Int) = weierstrass_model(hirzebruch_surface(NormalToricVariety, r); completeness_check = false)


@doc raw"""
    hypersurface_model_over_hirzebruch_surface(r::Int)

This method constructs a hypersurface model over a Hirzebruch surface.

# Examples
```jldoctest
julia> hypersurface_model_over_hirzebruch_surface(1)
Hypersurface model over a concrete base
```
"""
hypersurface_model_over_hirzebruch_surface(r::Int) = hypersurface_model(hirzebruch_surface(NormalToricVariety, r); completeness_check = false)


#####################################################
# 3. Models over del pezzo surface
#####################################################

@doc raw"""
    global_tate_model_over_del_pezzo_surface(b::Int)

This method constructs a global Tate model over a del-Pezzo surface.

# Examples
```jldoctest
julia> global_tate_model_over_del_pezzo_surface(3)
Global Tate model over a concrete base
```
"""
global_tate_model_over_del_pezzo_surface(b::Int) = global_tate_model(del_pezzo_surface(NormalToricVariety, b); completeness_check = false)


@doc raw"""
    weierstrass_model_over_del_pezzo_surface(b::Int)

This method constructs a Weierstrass model over a del-Pezzo surface.

# Examples
```jldoctest
julia> weierstrass_model_over_del_pezzo_surface(3)
Weierstrass model over a concrete base
```
"""
weierstrass_model_over_del_pezzo_surface(b::Int) = weierstrass_model(del_pezzo_surface(NormalToricVariety, b); completeness_check = false)


@doc raw"""
    hypersurface_model_over_del_pezzo_surface(b::Int)

This method constructs a hypersurface model over a del-Pezzo surface.

# Examples
```jldoctest
julia> hypersurface_model_over_del_pezzo_surface(3)
Hypersurface model over a concrete base
```
"""
hypersurface_model_over_del_pezzo_surface(b::Int) = hypersurface_model(del_pezzo_surface(NormalToricVariety, b); completeness_check = false)


#####################################################
# 4. Literature models
#####################################################

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
Multivariate polynomial ring in 10 variables over QQ graded by
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
Weierstrass model over a not fully specified base
```
"""
su5_weierstrass_model_over_arbitrary_3d_base() = weierstrass_model(su5_tate_model_over_arbitrary_3d_base())
