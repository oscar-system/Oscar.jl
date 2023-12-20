###################################################################
###################################################################
# 1: Attributes that work the same tor toric and non-toric settings
###################################################################
###################################################################


#####################################################
# 1.1 Hypersurface equation
#####################################################

@doc raw"""
    hypersurface_equation(h::HypersurfaceModel)

Return the hypersurface equation.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> hypersurface_equation(h);
```
"""
hypersurface_equation(h::HypersurfaceModel) = h.hypersurface_equation


#####################################################
# 1.2 Base, ambient space and fiber ambient space
#####################################################

@doc raw"""
    base_space(h::HypersurfaceModel)

Return the base space of the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> base_space(h)
Normal toric variety
```
"""
function base_space(h::HypersurfaceModel)
  base_fully_specified(h) || @vprint :HypersurfaceModel 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return h.base_space
end


@doc raw"""
    ambient_space(h::HypersurfaceModel)

Return the ambient space of the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> ambient_space(h)
Normal toric variety without torusfactor
```
"""
function ambient_space(h::HypersurfaceModel)
  base_fully_specified(h) || @vprint :HypersurfaceModel 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return h.ambient_space
end


@doc raw"""
    fiber_ambient_space(HypersurfaceModel)

Return the fiber ambient space of the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> fiber_ambient_space(h)
Normal toric variety
```
"""
fiber_ambient_space(h::HypersurfaceModel) = h.fiber_ambient_space



#########################################################################
#########################################################################
# 2: Attributes that define the corresponding Weierstrass and Tate models
#########################################################################
#########################################################################

@doc raw"""
    weierstrass_model(h::HypersurfaceModel)

Return the Weierstrass model corresponding to the
hypersurface model, provided that the latter is known.
"""
function weierstrass_model(h::HypersurfaceModel)
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  return get_attribute(h, :weierstrass_model)
end


@doc raw"""
    global_tate_model(h::HypersurfaceModel)

Return the global Tate model corresponding to the
hypersurface model, provided that the latter is known.
"""
function global_tate_model(h::HypersurfaceModel)
  @req has_attribute(h, :global_tate_model) "No corresponding global Tate model is known"
  return get_attribute(h, :global_tate_model)
end



############################################################################################################
############################################################################################################
# 3: Attributes that rest on the corresponding Weierstrass model and (currently) only work in toric settings
############################################################################################################
############################################################################################################


#####################################################
# 3.1 Calabi-Yau hypersurface
#####################################################

@doc raw"""
    calabi_yau_hypersurface(h::HypersurfaceModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> calabi_yau_hypersurface(h)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(h::HypersurfaceModel)
  @req typeof(base_space(h)) <: NormalToricVariety "Calabi-Yau hypersurface currently only supported for toric varieties/schemes as base space"
  base_fully_specified(h) || @vprint :HypersurfaceModel 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subvariety_of_toric_variety(ambient_space(h), [hypersurface_equation(h)])
end


#####################################################
# 3.2 Discriminant and singular loci
#####################################################

@doc raw"""
    discriminant(h::HypersurfaceModel)

Return the discriminant of the hypersurface model.
"""
@attr MPolyRingElem function discriminant(h::HypersurfaceModel)
  @req typeof(base_space(h)) <: NormalToricVariety "Discriminant currently only supported for toric varieties/schemes as base space"
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  return discriminant(weierstrass_model(h))
end


@doc raw"""
    singular_loci(h::HypersurfaceModel)

Return the singular loci of the hypersurface model, along with the order of
vanishing of the Weierstrass sections and discriminant ``(f, g, \Delta)```
at each locus. Also the refined Tate fiber type is returned.
"""
@attr Vector{<:Tuple{<:MPolyIdeal{<:MPolyRingElem}, Tuple{Int64, Int64, Int64}, String}} function singular_loci(h::HypersurfaceModel)
  @req typeof(base_space(h)) <: NormalToricVariety "Singular loci currently only supported for toric varieties/schemes as base space"
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  return singular_loci(weierstrass_model(h))
end
