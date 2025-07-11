##########################################
### (1) General attributes
##########################################

@doc raw"""
    base_space(m::AbstractFTheoryModel)

Return the base space of the F-theory model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> base_space(m)
Family of spaces of dimension d = 3
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
Family of spaces of dimension d = 5
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
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> fiber_ambient_space(t)
Normal toric variety
```
"""
function fiber_ambient_space(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :fiber_ambient_space) "fiber_ambient_space not supported for this F-theory model"
  return m.fiber_ambient_space
end




##########################################
### (2) Section attributes
##########################################

@doc raw"""
    defining_classes(m::AbstractFTheoryModel)

Returns the defining divisor classes of the given F-theory model.

For a detailed explanation of the concept, see [Literature Models](@ref).

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

Returns the model sections of the given F-theory model.

For a detailed explanation of the concept, see [Literature Models](@ref).

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

Returns the tunable sections of the given F-theory model.

For a detailed explanation of the concept, see [Literature Models](@ref).

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

Returns a dictionary that defines how the structural sections of the given F-theory model
are parametrized by the tunable sections.

For a detailed explanation of the concept, see [Literature Models](@ref).

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

Returns the divisor classes of the tunable sections expressed in the basis of the anticanonical divisor
of the base and the defining classes of the given model.

For a detailed explanation of the concept, see [Literature Models](@ref).

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

Returns the divisor classes of all model sections of the given F-theory model.

For a detailed explanation of the concept, see [Literature Models](@ref).

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

Returns the explicit polynomial expressions for all model sections of the given F-theory model.

For a detailed explanation of the concept, see [Literature Models](@ref).

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> explicit_model_sections(t)
Dict{String, QQMPolyRingElem} with 9 entries:
  "a6"  => 0
  "a21" => a21
  "a3"  => w^2*a32
  "w"   => w
  "a2"  => w*a21
  "a1"  => a1
  "a4"  => w^3*a43
  "a43" => a43
  "a32" => a32

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t2 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> explicit_model_sections(t2);
```
"""
function explicit_model_sections(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :explicit_model_sections) "explicit_model_sections not supported for this F-theory model"
  return m.explicit_model_sections
end




##############################################
### (3) Position of model in our database
##############################################

@doc raw"""
    model_index(m::AbstractFTheoryModel)
Return database index of a literature model. This index is a unique identifier that can be used to more conveniently construct the model. 
All models have a model_index and these will not change in the future.

```jldoctest
julia> t = literature_model(31)
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base -- F-theory weierstrass model dual to hypersurface model with fiber ambient space F_10 based on arXiv paper 1408.4808 Eq. (3.130)

julia> model_index(t)
31
```
"""
function model_index(m::AbstractFTheoryModel)
  directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  return parse(Int, model_indices["model" * literature_identifier(m) * ".json"])
end




##########################################
### (4) Data attributes macro
##########################################

macro define_model_attribute_getter(arg_expr)
    if !(arg_expr isa Expr && arg_expr.head == :tuple && length(arg_expr.args) == 2)
        error("Expected input like: (function_name, ReturnType)")
    end

    fname = arg_expr.args[1]
    rettype = arg_expr.args[2]
    sym = QuoteNode(fname)
    msg = "No $(replace(string(fname), '_' => ' ')) known for this model"

    doc = """
    Returns `$(fname)` of the F-theory model if known, otherwise throws an error.

    See [Literature Models](@ref) for more details.
    """

    return quote
        @doc $doc
        function $(esc(fname))(m::AbstractFTheoryModel)
            @req has_attribute(m, $sym) $msg
            return get_attribute(m, $sym)::$(esc(rettype))
        end
    end
end




##########################################
### (5) Meta data getters
##########################################

meta_data = [
    (:associated_literature_models, Vector{String}),
    (:arxiv_doi, String),
    (:arxiv_id, String),
    (:arxiv_link, String),
    (:arxiv_model_equation_number, String),
    (:arxiv_model_page, String),
    (:arxiv_model_section, String),
    (:arxiv_version, String),
    (:birational_literature_models, Vector{String}),
    (:journal_doi, String),
    (:journal_link, String),
    (:journal_model_equation_number, String),
    (:journal_model_page, String),
    (:journal_model_section, String),
    (:journal_pages, String),
    (:journal_report_numbers, Vector{String}),
    (:journal_volume, String),
    (:journal_name, String),
    (:journal_year, String),
    (:literature_identifier, String),
    (:model_description, String),
    (:model_parameters, Dict{String, Int64}),
    (:paper_authors, Vector{String}),
    (:paper_buzzwords, Vector{String}),
    (:paper_description, String),
    (:paper_title, String)
]

for (name, typ) in meta_data
    @eval @define_model_attribute_getter ($name, $typ)
end




##########################################
### (5) Predefined geometric data getters
##########################################

# Return a list of the known Mordell–Weil generating sections of the given model.  If no generating sections are known, an error is raised.
@define_model_attribute_getter (generating_sections, Vector{Vector{T}} where T <: Any)

# Return the list of all known resolutions for the given model. If no resolutions are known, an error is raised.
@define_model_attribute_getter (resolutions, Vector{Tuple{Vector{Vector{String}}, Vector{String}}})

# Return a list of lists of known Mordell–Weil generating sections for the given model after each known resolution. Each element of the outer list corresponds to a known resolution (in the same order), and each element of the list associated to a given resolution corresponds to a known generating section (in the same order). If no resolution generating sections are known, an error is raised.
@define_model_attribute_getter (resolution_generating_sections, Vector{Vector{Vector{Vector{T}}}} where T <: Any)

# Return a list of known Mordell–Weil zero sections for the given model after each known resolution. Each element of the list corresponds to a known resolution (in the same order). If no resolution zero sections are known, an error is raised.
@define_model_attribute_getter (resolution_zero_sections, Vector{Vector{Vector{T}}} where T <: Any)

# Return the list of all known weighted resolutions for the given model. If no weighted resolutions are known, an error is raised.
@define_model_attribute_getter (weighted_resolutions, Vector{Tuple{Vector{Tuple{Vector{String}, Vector{Int64}}}, Vector{String}}})

# Return a list of lists of known Mordell–Weil generating sections for the given model after each known weighted resolution. Each element of the outer list corresponds to a known weighted resolution (in the same order), and each element of the list associated to a given weighted resolution corresponds to a known generating section (in the same order). If no weighted resolution generating sections are known, an error is raised.
@define_model_attribute_getter (weighted_resolution_generating_sections, Vector{Vector{Vector{Vector{T}}}} where T <: Any)

# Return a list of known Mordell–Weil zero sections for the given model after each known weighted resolution. Each element of the list corresponds to a known weighted resolution (in the same order). If no weighted resolution zero sections are known, an error is raised.
@define_model_attribute_getter (weighted_resolution_zero_sections, Vector{Vector{Vector{T}}} where T <: Any)

# Return the zero section of the given model. If no zero section is known, an error is raised. This information is not typically stored as an attribute for Weierstrass and global Tate models, whose zero sections are known.
@define_model_attribute_getter (zero_section, Vector{T} where T <: Any)

# Return the zero section class of a model as a cohomology class in the toric ambient space. If no zero section class is known, an error is raised. This information is always available for Weierstrass and global Tate models, whose zero section classes are known.
@define_model_attribute_getter (zero_section_class, CohomologyClass)

# Return the index of the generator of the Cox ring of the ambient space, whose corresponding vanishing locus defines the zero section of a model. If no zero section class is known, an error is raised. This attribute is always set simultaneously with zero_section_class. This information is always available for Weierstrass and global Tate models, whose zero section classes are known.
@define_model_attribute_getter (zero_section_index, Int)

# Return the cohomology classes of the exceptional toric divisors of a model as a vector of cohomology classes in the toric ambient space. This information is only supported for models over a concrete base that is a normal toric variety, but is always available in this case. After a toric blow up this information is updated.
function exceptional_classes(m::AbstractFTheoryModel)
  @req base_space(m) isa NormalToricVariety "Exceptional divisor classes are only supported for models over a concrete base"
  return get_attribute(m, :exceptional_classes, Vector{CohomologyClass}())
end

# Return the indices of the generators of the Cox ring of the ambient space which correspond to exceptional divisors. This information is only supported for models over a concrete base that is a normal toric variety, but is always available in this case. After a toric blow up this information is updated.
@attr Vector{Int} function exceptional_divisor_indices(m::AbstractFTheoryModel)
  @req base_space(m) isa NormalToricVariety "Exceptional divisor indices are only supported for models over a concrete base"
  return get_attribute(m, :exceptional_divisor_indices, Vector{Int64}())
end

# Return the torsion sections of the given model. If no torsion sections are known, an error is raised.
@define_model_attribute_getter (torsion_sections, Vector{Vector{T}} where T <: Any)

# Return the gauge algebra of the given model. If no gauge algebra is known, an error is raised. This information is typically available for all models, however.
@define_model_attribute_getter (gauge_algebra, DirectSumLieAlgebra{QQBarFieldElem})

# Return list of lists of matrices, where each list of matrices corresponds to a gauge factor of the same index given by `gauge_algebra(m)`. These matrices are elements of the center of the corresponding gauge factor and quotienting by them replicates the action of some discrete group on the center of the lie algebra. This list combined with `gauge_algebra(m)` completely determines the gauge group of the model. If no gauge quotients are known, an error is raised.
@define_model_attribute_getter (global_gauge_group_quotient, Vector{Vector{String}})




##########################################
### (5) Computed geometric data getters
##########################################

@doc raw"""
    chern_class(m::AbstractFTheoryModel, k::Int; check::Bool = true)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
this very cohomology class ``h``, otherwise it raises an error.

The theory guarantees that the implemented algorithm works for toric ambient spaces which are
smooth and complete. The check for completeness can be very time consuming. This check can
be switched off by setting the optional argument `check` to the value `false`, as demonstrated below.

!!!warning
    This method works ONLY for F-theory models which are hypersurfaces in a toric ambient space.
  
!!!warning
    This method represents the Chern classes of said hypersurface by cohomology classes on the toric ambient space.
    These classes counterparts must be restricted to the hypersurface to truly represent the Chern class in question.
    Internally, we integrate those ambient space classes against the class of the hypersurface, which automatically
    executes the restriction to the hypersurface.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = chern_class(qsm_model, 4; check = false);

julia> is_trivial(h)
false
```
"""
function chern_class(m::AbstractFTheoryModel, k::Int; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"

  # Consistency checks
  @req k >= 0 "Chern class index must be non-negative"
  @req k <= dim(ambient_space(m)) - 1 "Chern class index must not exceed dimension of the space"

  # If thus far, no non-trivial Chern classes have been computed for this toric variety, add an "empty" vector
  if !has_attribute(m, :chern_classes)
    cs = Dict{Int64, CohomologyClass}()
    cs[0] = cohomology_class(ambient_space(m), one(cohomology_ring(ambient_space(m), check = check)), quick = true)
    diff = degree(leading_term(hypersurface_equation(m))) - sum(cox_ring(ambient_space(m)).d)
    cs[1] = cohomology_class(toric_divisor_class(ambient_space(m), diff), quick = true)
    set_attribute!(m, :chern_classes, cs)
    if k == 0
      return cs[0]
    elseif k == 1
      return cs[1]
    end
  end

  # Check if the Chern class in question is known
  cs = get_attribute(m, :chern_classes)::Dict{Int64, CohomologyClass}
  if haskey(cs, k)
    return cs[k]
  end

  # Check if we can compute the Chern classes for the toric ambient space
  if check
    @req is_smooth(ambient_space(m)) && is_complete(ambient_space(m)) "The Chern classes of the tangent bundle of the toric ambient space are only supported if the toric ambient space is smooth and complete"
  end

  # Chern class is not known, so compute and return it...
  cy = cohomology_class(toric_divisor_class(ambient_space(m), degree(leading_term(hypersurface_equation(m)))), quick = true)
  ck_ambient = chern_class(ambient_space(m), k, check = check)
  ckm1 = chern_class(m, k-1, check = check)
  new_poly = lift(polynomial(ck_ambient)) - lift(polynomial(cy)) * lift(polynomial(ckm1))
  coho_R = cohomology_ring(ambient_space(m), check = check)
  cs[k] = cohomology_class(ambient_space(m), coho_R(new_poly), quick = true)
  set_attribute!(m, :chern_classes, cs)
  return cs[k]
end


@doc raw"""
    chern_classes(m::AbstractFTheoryModel; check::Bool = true)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
a vector with the cohomology classes corresponding to all non-trivial Chern classes ``c_k`` of
``Y_n``. Otherwise, this methods raises an error.

As of right now, this method is computationally expensive for involved toric ambient spaces,
such as in the example below.

The theory guarantees that the implemented algorithm works for toric ambient spaces which are
simplicial and complete. The check for completeness can be very time consuming. This check can
be switched off by setting the optional argument `check` to the value `false`, as demonstrated below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = chern_classes(qsm_model; check = false);

julia> is_one(polynomial(h[0]))
true

julia> is_trivial(h[1])
true

julia> is_trivial(h[2])
false
```
"""
function chern_classes(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"
  for k in 0:dim(ambient_space(m))-1
    chern_class(m, k; check = check)
  end
  return get_attribute(m, :chern_classes)::Dict{Int64, CohomologyClass}
end


@doc raw"""
    euler_characteristic(m::AbstractFTheoryModel; check::Bool = true)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute the Euler characteristic. If this
assumptions is satisfied, this method returns the Euler characteristic, otherwise it raises an
error.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = euler_characteristic(qsm_model; check = false)
378

julia> h = euler_characteristic(qsm_model; check = false)
378
```
"""
@attr Int function euler_characteristic(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Euler characteristic of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Euler characteristic of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Euler characteristic of F-theory model currently supported only for toric ambient space"

  # Trigger potential short-cut computation of cohomology ring
  cohomology_ring(ambient_space(m); check)

  # Compute the cohomology class corresponding to the hypersurface equation
  cy = cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m))))

  # Compute the Euler characteristic
  return Int(integrate(chern_class(m, 4; check) * cy; check))
end




######################################################
### (6) Predefined geometric data getters for the QSMs
######################################################


### (4.1) Attributes regarding the polytope in the Kreuzer-Skarke database

@doc raw"""
    vertices(m::AbstractFTheoryModel)

This method returns the vertices of the polytope the the base of the F-theory QSM is
build from. Note that those vertices are normalized according to the Polymake standard
to rational numbers.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> vertices(qsm_model)
4-element Vector{Vector{QQFieldElem}}:
 [-1, -1, -1]
 [1, -1//2, -1//2]
 [-1, 2, -1]
 [-1, -1, 5]
```
"""
function vertices(m::AbstractFTheoryModel)
  @req has_attribute(m, :vertices) "No vertices known for this model"
  return get_attribute(m, :vertices)::Vector{Vector{QQFieldElem}}
end


@doc raw"""
    polytope_index(m::AbstractFTheoryModel)

Of the 3-dimensional reflexive polytope that the base of this F-theory model is build from,
this method returns the index within the Kreuzer-Skarke list.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> polytope_index(qsm_model)
4
```
"""
function polytope_index(m::AbstractFTheoryModel)
  @req has_attribute(m, :poly_index) "No polytope index known for this model"
  return get_attribute(m, :poly_index)::Int
end


@doc raw"""
    has_quick_triangulation(m::AbstractFTheoryModel)

For a 3-dimensional reflexive polytope in the Kreuzer-Skarke list, the list of
full (sometimes also called fine), regular, star triangulations can be extremely
large. Consequently, one may wonder if the triangulations can be enumerated in a
somewhat reasonable time (say 5 minutes on a personal computer). This method tries
to provide an answer to this. It returns `true` if one should expect a timely response
to the attempt to enumerate all (full, regular, star) triangulations. Otherwise, this
method returns `false`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> has_quick_triangulation(qsm_model)
true
```
"""
function has_quick_triangulation(m::AbstractFTheoryModel)
  @req has_attribute(m, :triang_quick) "It is not known if the base of this model can be triangulated quickly"
  return get_attribute(m, :triang_quick)::Bool
end


@doc raw"""
    max_lattice_pts_in_facet(m::AbstractFTheoryModel)

In order to enumerate the number of full, regular, star triangulations of a
3-dimensional reflexive polytope, it is possible to first find the corresponding
triangulations of all facets of the polytope [HT17](@cite). A first indication for
the complexity of this triangulation task is the maximum number of lattice points
in a facet of the polytope in question. This method returns this maximal number of
lattice points.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> max_lattice_pts_in_facet(qsm_model)
16
```
"""
function max_lattice_pts_in_facet(m::AbstractFTheoryModel)
  @req has_attribute(m, :max_lattice_pts_in_facet) "Maximal number of lattice points of facets not known for this model"
  return get_attribute(m, :max_lattice_pts_in_facet)::Int
end


@doc raw"""
    estimated_number_of_triangulations(m::AbstractFTheoryModel)

This method returns an estimate for the number of full, regular, star triangulations
of the 3-dimensional reflexive polytope, those triangulations define the possible base
spaces of the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> estimated_number_of_triangulations(qsm_model)
212533333333
```
"""
function estimated_number_of_triangulations(m::AbstractFTheoryModel)
  @req has_attribute(m, :estimated_number_of_triangulations) "Estimated number of (full, regular, star) triangulation not known for this model"
  return get_attribute(m, :estimated_number_of_triangulations)::Int
end


### (4.2) Attributes regarding the polytope in the Kreuzer-Skarke database


@doc raw"""
    kbar3(m::AbstractFTheoryModel)

Let Kbar denote the anticanonical class of the 3-dimensional base space of the F-theory QSM.
Of ample importance is the triple intersection number of Kbar, i.e. Kbar * Kbar * Kbar.
This method returns this intersection number.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> kbar3(qsm_model)
6
```
"""
function kbar3(m::AbstractFTheoryModel)
  @req has_attribute(m, :Kbar3) "Kbar3 not known for this model"
  return get_attribute(m, :Kbar3)::Int
end


@doc raw"""
    hodge_h11(m::AbstractFTheoryModel)

This methods return the Hodge number h11 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h11(qsm_model)
31
```
"""
function hodge_h11(m::AbstractFTheoryModel)
  @req has_attribute(m, :h11) "Hodge number h11 of ambient space not known for this model"
  return get_attribute(m, :h11)::Int
end


@doc raw"""
    hodge_h12(m::AbstractFTheoryModel)

This methods return the Hodge number h12 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h12(qsm_model)
10
```
"""
function hodge_h12(m::AbstractFTheoryModel)
  @req has_attribute(m, :h12) "Hodge number h12 of ambient space not known for this model"
  return get_attribute(m, :h12)::Int
end


@doc raw"""
    hodge_h13(m::AbstractFTheoryModel)

This methods return the Hodge number h13 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h13(qsm_model)
34
```
"""
function hodge_h13(m::AbstractFTheoryModel)
  @req has_attribute(m, :h13) "Hodge number h13 of ambient space not known for this model"
  return get_attribute(m, :h13)::Int
end


@doc raw"""
    hodge_h22(m::AbstractFTheoryModel)

This methods return the Hodge number h22 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h22(qsm_model)
284
```
"""
function hodge_h22(m::AbstractFTheoryModel)
  @req has_attribute(m, :h22) "Hodge number h22 of ambient space not known for this model"
  return get_attribute(m, :h22)::Int
end


### (4.3) Attributes regarding the Ci-curves


@doc raw"""
    genera_of_ci_curves(m::AbstractFTheoryModel)

This methods return the genera of the Ci curves.
Recall that Ci = V(xi, s), where xi is a homogeneous
coordinate of the 3-dimensional toric base space B3 of the
QSM hypersurface model in question, and s is a generic
section of the anticanonical bundle of B3. Consequently,
we may use the coordinates xi as labels for the curves Ci.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> keys_list = collect(keys(genera_of_ci_curves(qsm_model)));

julia> my_key = only(filter(k -> string(k) == "x7", keys_list))
x7

julia> genera_of_ci_curves(qsm_model)[my_key]
0
```
"""
function genera_of_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :genus_ci) "Genera of Ci curves not known for this model"
  return get_attribute(m, :genus_ci)
end


@doc raw"""
    degrees_of_kbar_restrictions_to_ci_curves(m::AbstractFTheoryModel)

The anticanonical divisor of the 3-dimensional toric base space B3 of the
QSM hypersurface model in question can be restricted to the Ci curves. The
result of this operation is a line bundle. This method returns the degree of
this line bundle for every Ci curve.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> keys_list = collect(keys(degrees_of_kbar_restrictions_to_ci_curves(qsm_model)));

julia> my_key = only(filter(k -> string(k) == "x7", keys_list))
x7

julia> degrees_of_kbar_restrictions_to_ci_curves(qsm_model)[my_key]
0
```
"""
function degrees_of_kbar_restrictions_to_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :degree_of_Kbar_of_tv_restricted_to_ci) "Degree of Kbar restriction to Ci curves not known for this model"
  return get_attribute(m, :degree_of_Kbar_of_tv_restricted_to_ci)
end


@doc raw"""
    topological_intersection_numbers_among_ci_curves(m::AbstractFTheoryModel)

The topological intersection numbers among Ci curves are also of ample importance.
This method returns those intersection numbers in the form of a matrix.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> n_rows(topological_intersection_numbers_among_ci_curves(qsm_model))
29

julia> n_columns(topological_intersection_numbers_among_ci_curves(qsm_model))
29
```
"""
function topological_intersection_numbers_among_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :intersection_number_among_ci_cj) "Topological intersection numbers among Ci curves not known for this model"
  return get_attribute(m, :intersection_number_among_ci_cj)
end


@doc raw"""
    indices_of_trivial_ci_curves(m::AbstractFTheoryModel)

Some of the Ci curves are trivial, in that V(xi, s) is the empty set.
This method returns the vector of all indices of trivial Ci curves.
That is, should V(x23, s) be the empty set, then 23 will be included in
the returned list.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> indices_of_trivial_ci_curves(qsm_model)
10-element Vector{Int64}:
 23
 22
 18
 19
 20
 26
 10
 11
 12
 15
```
"""
function indices_of_trivial_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :index_facet_interior_divisors) "Degree of Kbar restriction to Ci curves not known for this model"
  return get_attribute(m, :index_facet_interior_divisors)::Vector{Int}
end


@doc raw"""
    topological_intersection_numbers_among_nontrivial_ci_curves(m::AbstractFTheoryModel)

The topological intersection numbers among the non-trivial Ci curves are used
frequently. This method returns those intersection numbers in the form of a matrix.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> n_rows(topological_intersection_numbers_among_nontrivial_ci_curves(qsm_model))
19

julia> n_columns(topological_intersection_numbers_among_nontrivial_ci_curves(qsm_model))
19
```
"""
function topological_intersection_numbers_among_nontrivial_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :intersection_number_among_nontrivial_ci_cj) "Topological intersection numbers among non-trivial Ci curves not known for this model"
  return get_attribute(m, :intersection_number_among_nontrivial_ci_cj)
end


### (4.4) Attributes regarding the dual graph


@doc raw"""
    dual_graph(m::AbstractFTheoryModel)

This method returns the dual graph of the QSM model in question.
Note that no labels are (currently) attached to the vertices/nodes or edges.
To understand/read this graph correctly, please use the methods listed below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> dual_graph(qsm_model)
Undirected graph with 21 nodes and the following edges:
(5, 1)(6, 5)(7, 6)(8, 7)(9, 4)(9, 8)(10, 1)(11, 4)(12, 3)(12, 10)(13, 3)(13, 11)(14, 1)(15, 4)(16, 3)(17, 3)(18, 2)(18, 14)(19, 2)(19, 15)(20, 2)(20, 16)(21, 2)(21, 17)
```
"""
function dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :dual_graph) "Dual graph not known for this model"
  return get_attribute(m, :dual_graph)
end


@doc raw"""
    components_of_dual_graph(m::AbstractFTheoryModel)

This method returns a vector with labels for each node/vertex of the dual graph of the QSM
model in question. Those labels allow to understand the geometric origin of the node/vertex.

Specifically, recall that those nodes are associated to the Ci-curves, which are in turn
given by Ci = V(xi, s). xi is a homogeneous coordinate of the 3-dimensional toric base space
B3 of the QSM in question, and s is a generic section of the anticanonical bundle of B3.

Only non-trivial Ci = V(xi, s) correspond to vertices/nodes of the dual graph.

If Ci = V(xi, s) is irreducible and corresponds to the k-th component, then the label "Ci"
appears at position k of the vector returned by this method. However, if Ci = V(xi, s) is
reducible, then we introduce the labels Ci-0, Ci-1, Ci-2 etc. for those irreducible
components of Ci. If Ci-0 corresponds to the k-th components of the dual graph,
then the label "Ci-0" appears at position k of the vector returned by this method.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> components_of_dual_graph(qsm_model)
21-element Vector{String}:
 "C0"
 "C1"
 "C2"
 "C3"
 "C4"
 "C5"
 "C6"
 "C7"
 "C8"
 "C9"
 ⋮
 "C16"
 "C17"
 "C21"
 "C24-0"
 "C24-1"
 "C25"
 "C27"
 "C28-0"
 "C28-1"
```
"""
function components_of_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :components_of_dual_graph) "Components of dual graph not known for this model"
  return get_attribute(m, :components_of_dual_graph)
end


@doc raw"""
    degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)

The anticanonical bundle of the toric 3-dimensional base space of the F-theory QSM in
question can be restricted to the (geometric counterparts of the) nodes/vertices of
the dual graph. The result is a line bundle for each node/vertex. This method returns
a vector with the degrees of these line bundles.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> degrees_of_kbar_restrictions_to_components_of_dual_graph(qsm_model)["C28-1"]
0
```
"""
function degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph) "Degree of Kbar restricted to components of dual graph not known for this model"
  return get_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph)
end


@doc raw"""
    genera_of_components_of_dual_graph(m::AbstractFTheoryModel)

This methods returns a vector with the genera of the (geometric
counterparts of the) nodes/vertices of the dual graph.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> genera_of_components_of_dual_graph(qsm_model)["C28-1"]
0
```
"""
function genera_of_components_of_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :genus_of_components_of_dual_graph) "Genera of components of dual graph not known for this model"
  return get_attribute(m, :genus_of_components_of_dual_graph)
end


### (4.5) Attributes regarding the simplified dual graph


@doc raw"""
    simplified_dual_graph(m::AbstractFTheoryModel)

This method returns the simplified dual graph of the QSM model in question.
Note that no labels are (currently) attached to the vertices/nodes or edges.
To understand/read this graph correctly, please use the methods listed below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> simplified_dual_graph(qsm_model)
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)
```
"""
function simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :simplified_dual_graph) "Simplified dual graph not known for this model"
  return get_attribute(m, :simplified_dual_graph)
end


@doc raw"""
    components_of_simplified_dual_graph(m::AbstractFTheoryModel)

This method returns a vector with labels for each node/vertex of the simplified dual graph.
Otherwise, works identical to `components_of_dual_graph(m::AbstractFTheoryModel)`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> components_of_simplified_dual_graph(qsm_model)
4-element Vector{String}:
 "C0"
 "C1"
 "C2"
 "C3"
```
"""
function components_of_simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :components_of_simplified_dual_graph) "Components of simplified dual graph not known for this model"
  return get_attribute(m, :components_of_simplified_dual_graph)
end


@doc raw"""
    degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(m::AbstractFTheoryModel)

Same as `degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)`,
but for the simplified dual graph.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(qsm_model)["C2"]
2
```
"""
function degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph) "Degree of Kbar restricted to components of simplified dual graph not known for this model"
  return get_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph)
end


@doc raw"""
    genera_of_components_of_simplified_dual_graph(m::AbstractFTheoryModel)

This methods returns a vector with the genera of the (geometric
counterparts of the) nodes/vertices of the dual graph.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> genera_of_components_of_simplified_dual_graph(qsm_model)["C2"]
0
```
"""
function genera_of_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :genus_of_components_of_simplified_dual_graph) "Genera of components of simplified dual graph not known for this model"
  return get_attribute(m, :genus_of_components_of_simplified_dual_graph)
end


######################################################################################
### (5) Attributes for flux families (not exported, rather for serialization overhaul)
######################################################################################

@attr QQMatrix function matrix_integral_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_integral(special_flux_family(m, check = check))
end

@attr QQMatrix function matrix_rational_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_rational(special_flux_family(m, check = check))
end

@attr Vector{QQFieldElem} function offset_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return offset(special_flux_family(m, check = check))
end

@attr QQMatrix function matrix_integral_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_integral(special_flux_family(m, not_breaking = true; check = check))
end

@attr QQMatrix function matrix_rational_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_rational(special_flux_family(m, not_breaking = true; check = check))
end

@attr Vector{QQFieldElem} function offset_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return offset(special_flux_family(m, not_breaking = true; check = check))
end
