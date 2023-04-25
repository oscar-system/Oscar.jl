#####################################################
# 1: Fiber analysis
#####################################################

@doc raw"""
    analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})

Determine the fiber of a (singular) global Tate model over a particular base locus.
```
"""
function analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
  
  # This method only works if the model is defined over a toric variety over toric scheme
  @req typeof(base_space(model)) <: ToricCoveredScheme "Analysis of fibers currently only supported for toric scheme/variety as base space"
  
  # Ideal of the defining polynomial
  hypersurface_ideal = ideal([tate_polynomial(model)])
  
  # Toric ambient space
  tas = ambient_space(model)
  
  # Various important ideals
  irr = irrelevant_ideal(tas);
  sri = stanley_reisner_ideal(tas);
  lin = ideal_of_linear_relations(tas);
  
  # Singular loci
  sing_loc = singular_loci(model)
  
  # Pick out the singular loci that are more singular than an I_1
  # Then keep only the locus and not the extra info about it
  interesting_singular_loci = map(tup -> tup[1], filter(locus -> locus[2][3] > 1, sing_loc))
  
  # This is a kludge to map polynomials on the base into the ambient space, and should be fixed once the ambient space constructors supply such a map
  base_coords = parent(gens(interesting_singular_loci[1])[1])
  ambient_coords = parent(tate_polynomial(model))
  base_to_ambient_ring_map = hom(base_coords, ambient_coords, gens(ambient_coords)[1:end-3])
  
  # Resolved model
  strict_transform, exceptionals, crepant, res_irr, res_sri, res_lin, res_S, res_S_gens, res_ring_map = _blowup_global_sequence(hypersurface_ideal, centers, irr, sri, lin)
  if !crepant
      @warn "The given sequence of blowups is not crepant"
  end
  
  loci_fiber_intersections = Tuple{MPolyIdeal{QQMPolyRingElem}, Vector{Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}}}[]
  for locus in interesting_singular_loci
    # Currently have to get the ungraded ideal generators by hand using .f
    ungraded_locus = ideal(map(gen -> base_to_ambient_ring_map(gen).f, gens(locus)))
    
    # Potential components of the fiber over this locus
    # For now, we only consider the associated prime ideal,
    # but we may later want to actually consider the primary ideals
    potential_components = map(pair -> pair[2], primary_decomposition(strict_transform + res_ring_map(ungraded_locus)))
    
    # Filter out the trivial loci among the potential components
    components = filter(component -> _is_nontrivial(component, res_irr), potential_components)
    
    # Check the pairwise intersections of the components
    intersections = Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}[]
    for i in 1:length(components) - 1
      for j in i + 1:length(components)
        intersection = filter(candidate_locus -> _is_nontrivial(candidate_locus, res_irr), map(pair -> pair[2], primary_decomposition(components[i] + components[j])))
        push!(intersections, ((i, j), intersection))
      end
    end
    
    push!(loci_fiber_intersections, (ungraded_locus, intersections))
  end
  
  return loci_fiber_intersections

end


#####################################################
# 2: Adjust the description for the global Tate model
#####################################################

@doc raw"""
    set_description(t::GlobalTateModel, description::String)

Set a description for a global Tate model.

```jldoctest
julia> t = literature_tate_model("1109.3454", "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> set_description(t, "An SU(5)xU(1) GUT-model")

julia> t
Global Tate model over a not fully specified base -- An SU(5)xU(1) GUT-model based on arxiv paper 1109.3454 (equ. 3.5)
```
"""
function set_description(t::GlobalTateModel, description::String)
  set_attribute!(t, :description => description)
end


#####################################################
# 3: Add a resolution
#####################################################

@doc raw"""
    add_resolution(t::GlobalTateModel, resolution::Vector{Vector{String}})

Set a description for a global Tate model.

```jldoctest
julia> t = literature_tate_model("1109.3454", "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> add_resolution(t, [["x", "y"], ["y", "s", "w"], ["s", "e4"], ["s", "e3"], ["s", "e1"], ["s", "w", "e3", "e1", "e2"]])

julia> length(resolutions(t))
2
```
"""
function add_resolution(t::GlobalTateModel, resolution::Vector{Vector{String}})
  if has_attribute(t, :resolutions)
    known_resolutions = resolutions(t)
    if (resolution in known_resolutions) == false
      push!(known_resolutions, resolution)
      set_attribute!(t, :resolutions => known_resolutions)
    end
  else
    set_attribute!(t, :resolutions => [resolution])
  end
end


#####################################################
# 4: Resolve a model with a known resolution
#####################################################

@doc raw"""
    resolve(t::GlobalTateModel, index::Int)

Resolve a global Tate model with the index-th resolution that is known.

Careful: Currently, this assumes that all blowups are toric blowups.
We hope to remove this requirement in the near future.

```jldoctest
julia> t = literature_tate_model("1109.3454", "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> v = resolve(t, 1)
Normal toric variety

julia> cox_ring(v)
Multivariate Polynomial Ring in 14 variables a10, a21, a32, a43, ..., s over Rational Field graded by 
  a10 -> [0 0 0 0 0 0]
  a21 -> [0 0 0 0 0 0]
  a32 -> [0 0 0 0 0 0]
  a43 -> [0 0 0 0 0 0]
  a65 -> [0 0 0 0 0 0]
  w -> [1 0 0 0 0 0]
  x -> [0 1 0 0 0 0]
  y -> [0 0 1 0 0 0]
  z -> [0 0 0 1 0 0]
  e1 -> [0 0 0 0 1 0]
  e4 -> [0 0 0 0 0 1]
  e2 -> [-1 -1 1 -1 -1 0]
  e3 -> [0 1 -1 1 0 -1]
  s -> [2 -1 0 2 1 1]
```
"""
function resolve(t::GlobalTateModel, index::Int)
  @req has_attribute(t, :resolutions) "No resolutions known for this model"
  @req index > 0 "The resolution must be specified by a non-negative integer"
  @req index < length(resolutions(t)) "The resolution must be specified by an integer that is not larger than the number of known resolutions"
  
  # Gather information for resolution
  resolution = resolutions(t)[index]
  nr_blowups = length(resolution)-1
  
  # Is this a sequence of toric blowups? (To be extended with @HechtiDerLachs and ToricSchemes).
  function string_to_poly(R, s::String)
    evaluate(eval(Meta.parse("_, ($(join(symbols(R), ','))) = PolynomialRing(ZZ, $(ngens(R)), cached = false);"*s)), gens(R));
  end
  resolved_ambient_space = toric_ambient_space(t)
  R, gR = PolynomialRing(QQ, vcat([string(g) for g in gens(cox_ring(resolved_ambient_space))], resolution[nr_blowups+1]))
  for k in 1:nr_blowups
    @req all(x -> x in gR, [string_to_poly(R, p) for p in resolution[k]]) "Blowup currently not supported"
  end
  
  # Perform resolution
  for k in 1:nr_blowups
    S = cox_ring(resolved_ambient_space)
    resolved_ambient_space = blow_up(resolved_ambient_space, ideal([string_to_poly(S, g) for g in resolution[k]]); coordinate_name = resolution[nr_blowups + 1][k], set_attributes = true)
  end
  return resolved_ambient_space
end
