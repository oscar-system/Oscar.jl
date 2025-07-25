@doc raw"""
    defining_classes(m::AbstractFTheoryModel)

Return the defining divisor classes of the given F-theory model.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> defining_classes(t)
Dict{String, ToricDivisorClass} with 1 entry:
  "w" => Divisor class on a normal toric variety
```
"""
function defining_classes(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :defining_classes) "defining_classes not supported for this F-theory model"
  return m.defining_classes
end


@doc raw"""
    model_sections(m::AbstractFTheoryModel)

Return the model sections of the given F-theory model.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> model_sections(m)
9-element Vector{String}:
 "a6"
 "a21"
 "a3"
 "w"
 "a2"
 "a1"
 "a4"
 "a43"
 "a32"
```
"""
model_sections(m::AbstractFTheoryModel) = collect(keys(explicit_model_sections(m)))


@doc raw"""
    tunable_sections(m::AbstractFTheoryModel)

Return the tunable sections of the given F-theory model.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tunable_sections(m)
5-element Vector{String}:
 "a21"
 "w"
 "a1"
 "a43"
 "a32"
```
"""
@attr Vector{String} tunable_sections(m::AbstractFTheoryModel) = collect(setdiff(keys(explicit_model_sections(m)), keys(model_section_parametrization(m))))


@doc raw"""
    model_section_parametrization(m::AbstractFTheoryModel)

Return a dictionary that defines how the structural sections of the given F-theory model
are parametrized by the tunable sections.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> model_section_parametrization(t)
Dict{String, MPolyRingElem} with 4 entries:
  "a6" => 0
  "a3" => w^2*a32
  "a2" => w*a21
  "a4" => w^3*a43
```
"""
function model_section_parametrization(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :model_section_parametrization) "model_section_parametrization not supported for this F-theory model"
  return m.model_section_parametrization
end


@doc raw"""
    classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(m::AbstractFTheoryModel)

Return the divisor classes of the tunable sections expressed in the basis of the anticanonical divisor
of the base and the defining classes of the given model.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(m)
Dict{String, Vector{Int64}} with 5 entries:
  "a21" => [2, -1]
  "w"   => [0, 1]
  "a1"  => [1, 0]
  "a43" => [4, -3]
  "a32" => [3, -2]
```
"""
classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(m::AbstractFTheoryModel) = get_attribute(m, :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes, "Classes of tunable sections in basis of Kbar and defining classes not known for this model")::Dict{String, Vector{Int64}}


@doc raw"""
    classes_of_model_sections(m::AbstractFTheoryModel)

Return the divisor classes of all model sections of the given F-theory model.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> classes_of_model_sections(t)
Dict{String, ToricDivisorClass} with 9 entries:
  "a21" => Divisor class on a normal toric variety
  "a6"  => Divisor class on a normal toric variety
  "a3"  => Divisor class on a normal toric variety
  "w"   => Divisor class on a normal toric variety
  "a2"  => Divisor class on a normal toric variety
  "a1"  => Divisor class on a normal toric variety
  "a43" => Divisor class on a normal toric variety
  "a4"  => Divisor class on a normal toric variety
  "a32" => Divisor class on a normal toric variety
```
"""
@attr Dict{String, ToricDivisorClass} function classes_of_model_sections(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :explicit_model_sections) "explicit_model_sections not supported for this F-theory model"
  @req base_space(m) isa NormalToricVariety "classes_of_model_sections only supported for models over a toric base"
  my_dict = Dict{String, ToricDivisorClass}()
  for (key, value) in explicit_model_sections(m)
    sec = value
    @req is_homogeneous(sec) "Encountered a non-homogeneous model section"
    if is_zero(sec)
      my_dict[key] = trivial_divisor_class(base_space(m))
    else
      deg = collect(keys(homogeneous_components(sec)))[1]
      my_dict[key] = toric_divisor_class(base_space(m), deg)
    end
  end
  return my_dict
end


@doc raw"""
    explicit_model_sections(m::AbstractFTheoryModel)

Return the explicit polynomial expressions for all model sections of the given F-theory model.

For a detailed explanation of the concept, see [Model Parameters, Defining Classes, and Related Concepts](@ref model_section_explanation_section).

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> explicit_model_sections(t)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 9 entries:
  "a6"  => 0
  "a21" => a21
  "a3"  => w^2*a32
  "w"   => w
  "a2"  => w*a21
  "a1"  => a1
  "a4"  => w^3*a43
  "a43" => a43
  "a32" => a32
```
"""
function explicit_model_sections(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :explicit_model_sections) "explicit_model_sections not supported for this F-theory model"
  return m.explicit_model_sections
end
