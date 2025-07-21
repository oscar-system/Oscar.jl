### (1) Attributes regarding the polytope in the Kreuzer-Skarke database

@doc raw"""
    vertices(m::AbstractFTheoryModel)

This method returns the vertices of the polytope the the base of the F-theory QSM is
build from. Note that those vertices are normalized according to the Polymake standard
to rational numbers.

# Examples
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

# Examples
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

# Examples
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

# Examples
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

# Examples
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


### (2) Attributes regarding the polytope in the Kreuzer-Skarke database


@doc raw"""
    kbar3(m::AbstractFTheoryModel)

Let Kbar denote the anticanonical class of the 3-dimensional base space of the F-theory QSM.
Of ample importance is the triple intersection number of Kbar, i.e. Kbar * Kbar * Kbar.
This method returns this intersection number.

# Examples
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

# Examples
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

# Examples
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

# Examples
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

# Examples
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


### (3) Attributes regarding the Ci-curves


@doc raw"""
    genera_of_ci_curves(m::AbstractFTheoryModel)

This methods return the genera of the Ci curves.
Recall that Ci = V(xi, s), where xi is a homogeneous
coordinate of the 3-dimensional toric base space B3 of the
QSM hypersurface model in question, and s is a generic
section of the anticanonical bundle of B3. Consequently,
we may use the coordinates xi as labels for the curves Ci.

# Examples
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

# Examples
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

# Examples
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

# Examples
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

# Examples
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


### (4) Attributes regarding the dual graph


@doc raw"""
    dual_graph(m::AbstractFTheoryModel)

This method returns the dual graph of the QSM model in question.
Note that no labels are (currently) attached to the vertices/nodes or edges.
To understand/read this graph correctly, please use the methods listed below.

# Examples
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

# Examples
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

# Examples
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

# Examples
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


### (5) Attributes regarding the simplified dual graph


@doc raw"""
    simplified_dual_graph(m::AbstractFTheoryModel)

This method returns the simplified dual graph of the QSM model in question.
Note that no labels are (currently) attached to the vertices/nodes or edges.
To understand/read this graph correctly, please use the methods listed below.

# Examples
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

# Examples
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

# Examples
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

# Examples
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
