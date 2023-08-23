#####################################################
# 1: Adjust the description for the model
#####################################################

@doc raw"""
    set_description(m::AbstractFTheoryModel, description::String)

Set a description for a model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> set_description(m, "An SU(5)xU(1) GUT-model")

julia> m
Global Tate model over a not fully specified base -- An SU(5)xU(1) GUT-model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function set_description(m::AbstractFTheoryModel, description::String)
  set_attribute!(m, :model_description => description)
end


#####################################################
# 2: Add a resolution
#####################################################

@doc raw"""
    add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})

Add a known resolution for a model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_resolution(m, [["x", "y"], ["y", "s", "w"], ["s", "e4"], ["s", "e3"], ["s", "e1"]], ["s", "w", "e3", "e1", "e2"])

julia> length(resolutions(m))
2
```
"""
function add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
  @req length(exceptionals) == length(centers) "Number of exceptionals must match number of centers"

  resolution = [centers, exceptionals]
  if has_attribute(m, :resolutions)
    known_resolutions = resolutions(m)
    if (resolution in known_resolutions) == false
      push!(known_resolutions, resolution)
      set_attribute!(m, :resolutions => known_resolutions)
    end
  else
    set_attribute!(m, :resolutions => [resolution])
  end
end


#####################################################
# 3: Resolve a model with a known resolution
#####################################################

@doc raw"""
    resolve(m::AbstractFTheoryModel, index::Int)

Resolve a model with the index-th resolution that is known.

Careful: Currently, this assumes that all blowups are toric blowups.
We hope to remove this requirement in the near future.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> v = resolve(m, 1)
Scheme of a toric variety

julia> cox_ring(v)
Multivariate polynomial ring in 13 variables over QQ graded by 
  a1 -> [1 0 0 0 0 0 0 0]
  a21 -> [0 1 0 0 0 0 0 0]
  a32 -> [-1 2 0 0 0 0 0 0]
  a43 -> [-2 3 0 0 0 0 0 0]
  w -> [0 0 1 0 0 0 0 0]
  x -> [0 0 0 1 0 0 0 0]
  y -> [0 0 0 0 1 0 0 0]
  z -> [0 0 0 0 0 1 0 0]
  e1 -> [0 0 0 0 0 0 1 0]
  e4 -> [0 0 0 0 0 0 0 1]
  e2 -> [1 -1 -1 -1 1 -1 -1 0]
  e3 -> [1 0 0 1 -1 1 0 -1]
  s -> [-2 2 2 -1 0 2 1 1]
```
"""
function resolve(m::AbstractFTheoryModel, index::Int)
  @req has_attribute(m, :resolutions) "No resolutions known for this model"
  @req index > 0 "The resolution must be specified by a non-negative integer"
  @req index <= length(resolutions(m)) "The resolution must be specified by an integer that is not larger than the number of known resolutions"
  
  # Gather information for resolution
  centers, exceptionals = resolutions(m)[index]
  nr_blowups = length(centers)
  
  # Is this a sequence of toric blowups? (To be extended with @HechtiDerLachs and ToricSchemes).
  resolved_ambient_space = underlying_toric_variety(ambient_space(m))
  R, gR = PolynomialRing(QQ, vcat([string(g) for g in gens(cox_ring(resolved_ambient_space))], exceptionals))
  for center in centers
    @req all(x -> x in gR, [eval_poly(p, R) for p in center]) "Non-toric blowup currently not supported"
  end
  
  # Perform resolution
  for k in 1:nr_blowups
    S = cox_ring(resolved_ambient_space)
    resolved_ambient_space = domain(blow_up(resolved_ambient_space, ideal([eval_poly(g, S) for g in centers[k]]); coordinate_name = exceptionals[k], set_attributes = true))
  end
  return toric_covered_scheme(resolved_ambient_space)
end
