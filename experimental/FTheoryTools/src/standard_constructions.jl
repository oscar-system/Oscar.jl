#####################################################
# 1. Models over projective space
#####################################################

@doc raw"""
    global_tate_model_over_projective_space(d::Int)

Construct a global Tate model over the ``d``-dimensional projective space,
represented as a toric variety. The Tate sections ``a_i`` are
automatically generated with pseudorandom coefficients.

# Examples
```jldoctest
julia> global_tate_model_over_projective_space(3)
Global Tate model over a concrete base
```
"""
global_tate_model_over_projective_space(d::Int) = global_tate_model(projective_space(NormalToricVariety, d); completeness_check = false)


@doc raw"""
    weierstrass_model_over_projective_space(d::Int)

Construct a Weierstrass model over the ``d``-dimensional projective space,
represented as a toric variety. The Weierstrass sections ``f`` and ``g`` are
automatically generated with pseudorandom coefficients.

# Examples
```jldoctest
julia> weierstrass_model_over_projective_space(3)
Weierstrass model over a concrete base
```
"""
weierstrass_model_over_projective_space(d::Int) = weierstrass_model(projective_space(NormalToricVariety, d); completeness_check = false)


#####################################################
# 2. Models over hirzebruch surfaces
#####################################################

@doc raw"""
    global_tate_model_over_hirzebruch_surface(r::Int)

Construct a global Tate model over the Hirzebruch surface ``F_r``,
represented as a toric variety. The Tate sections ``a_i`` are
automatically generated with pseudorandom coefficients.

# Examples
```jldoctest
julia> global_tate_model_over_hirzebruch_surface(1)
Global Tate model over a concrete base
```
"""
global_tate_model_over_hirzebruch_surface(r::Int) = global_tate_model(hirzebruch_surface(NormalToricVariety, r); completeness_check = false)


@doc raw"""
    weierstrass_model_over_hirzebruch_surface(r::Int)

Construct a Weierstrass model over the Hirzebruch surface ``F_r``,
represented as a toric variety. The Weierstrass sections ``f`` and ``g`` are
automatically generated with pseudorandom coefficients.

# Examples
```jldoctest
julia> weierstrass_model_over_hirzebruch_surface(1)
Weierstrass model over a concrete base
```
"""
weierstrass_model_over_hirzebruch_surface(r::Int) = weierstrass_model(hirzebruch_surface(NormalToricVariety, r); completeness_check = false)


#####################################################
# 3. Models over del pezzo surface
#####################################################

@doc raw"""
    global_tate_model_over_del_pezzo_surface(b::Int)

Construct a global Tate model over the del Pezzo surface ``\text{dP}_b``,
represented as a toric variety. The Tate sections ``a_i`` are
automatically generated with pseudorandom coefficients.

# Examples
```jldoctest
julia> global_tate_model_over_del_pezzo_surface(3)
Global Tate model over a concrete base
```
"""
global_tate_model_over_del_pezzo_surface(b::Int) = global_tate_model(del_pezzo_surface(NormalToricVariety, b); completeness_check = false)


@doc raw"""
    weierstrass_model_over_del_pezzo_surface(b::Int)

Construct a Weierstrass model over the del Pezzo surface ``\text{dP}_b``,
represented as a toric variety. The Weierstrass sections ``f`` and ``g`` are
automatically generated with pseudorandom coefficients.

# Examples
```jldoctest
julia> weierstrass_model_over_del_pezzo_surface(3)
Weierstrass model over a concrete base
```
"""
weierstrass_model_over_del_pezzo_surface(b::Int) = weierstrass_model(del_pezzo_surface(NormalToricVariety, b); completeness_check = false)
