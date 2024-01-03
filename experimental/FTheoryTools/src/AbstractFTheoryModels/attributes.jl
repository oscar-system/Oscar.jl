@doc raw"""
    base_space(m::AbstractFTheoryModel)

Return the base space of the F-theory model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> base_space(m)
A family of spaces of dimension d = 3
```
"""
function base_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return m.base_space
end


@doc raw"""
    ambient_space(m::AbstractFTheoryModel)

Return the ambient space of the F-theory model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> ambient_space(m)
A family of spaces of dimension d = 5
```
"""
function ambient_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return m.ambient_space
end


@doc raw"""
    fiber_ambient_space(m::AbstractFTheoryModel)

Return the fiber ambient space of an F-theory model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base

julia> fiber_ambient_space(t)
Normal toric variety
```
"""
function fiber_ambient_space(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :fiber_ambient_space) "fiber_ambient_space not supported for this F-theory model"
  return m.fiber_ambient_space
end
