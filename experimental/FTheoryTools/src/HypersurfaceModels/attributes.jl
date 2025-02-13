###################################################################
###################################################################
# 1: Attributes that work the same tor toric and non-toric settings
###################################################################
###################################################################

@doc raw"""
    hypersurface_equation(h::HypersurfaceModel)

Return the hypersurface equation.

```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, defining_classes = Dict("b" => b))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> hypersurface_equation(h);
```
"""
hypersurface_equation(h::HypersurfaceModel) = h.hypersurface_equation


@doc raw"""
    hypersurface_equation_parametrization(h::HypersurfaceModel)

Return the parametrization of the hypersurface
equation by the model sections.

```jldoctest
julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> explicit_model_sections(h)
Dict{String, MPolyRingElem} with 5 entries:
  "c2" => c2
  "c1" => c1
  "c3" => c3
  "b"  => b
  "c0" => c0

julia> hypersurface_equation_parametrization(h)
b*w*v^2 - c0*u^4 - c1*u^3*v - c2*u^2*v^2 - c3*u*v^3 + w^2
```
"""
hypersurface_equation_parametrization(h::HypersurfaceModel) = h.hypersurface_equation_parametrization



#########################################################################
#########################################################################
# 2: Attributes that define the corresponding Weierstrass and Tate models
#########################################################################
#########################################################################

@doc raw"""
    weierstrass_model(h::HypersurfaceModel)

Return the Weierstrass model corresponding to the
hypersurface model, provided that the former is known.

```jldoctest
julia> t = literature_model(14)
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> weierstrass_model(t)
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base -- F-theory weierstrass model dual to hypersurface model with fiber ambient space F_1 based on arXiv paper 1408.4808 Eq. (3.4)
```
"""
function weierstrass_model(h::HypersurfaceModel)
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  w = get_attribute(h, :weierstrass_model)
  if w isa String
    directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
    model_indices = JSON.parsefile(directory * "model_indices.json")
    if is_base_space_fully_specified(h)
      w_model = literature_model(parse(Int, model_indices[w]), base_space = base_space(h), defining_classes = defining_classes(h), completeness_check = false)
    else
      w_model = literature_model(parse(Int, model_indices[w]))
    end
    set_weierstrass_model(h, w_model)
    return w_model
  else
    return w
  end
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
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, defining_classes = Dict("b" => b))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> calabi_yau_hypersurface(h)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(h::HypersurfaceModel)
  @req base_space(h) isa NormalToricVariety "Calabi-Yau hypersurface currently only supported for toric varieties as base space"
  is_base_space_fully_specified(h) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
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
  @req base_space(h) isa NormalToricVariety "Discriminant currently only supported for toric varieties as base space"
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
  @req base_space(h) isa NormalToricVariety "Singular loci currently only supported for toric varieties as base space"
  @req has_attribute(h, :weierstrass_model) || has_attribute(h, :global_tate_model) "No corresponding Weierstrass model or global Tate model is known"
  return has_attribute(h, :weierstrass_model) ? singular_loci(weierstrass_model(h)) : singular_loci(global_tate_model(h))
end
