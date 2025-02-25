#####################################################
# Setters
#####################################################

@doc raw"""
    set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)

Allows to define the Weierstrass model corresponding to the hypersurface model.
"""
function set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)
  set_attribute!(h, :weierstrass_model => w)
end

@doc raw"""
    set_global_tate_model(h::HypersurfaceModel, w::GlobalTateModel)

Allows to define the global Tate model corresponding to the hypersurface model.
"""
function set_global_tate_model(h::HypersurfaceModel, t::GlobalTateModel)
  set_attribute!(h, :global_tate_model => t)
end



#####################################################
# 2: Tune a Hypersurface model
#####################################################

@doc raw"""
    function tune(h::HypersurfaceModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)

Tune a hypersurface model by fixing a special choice for the model sections.
Note that it is in particular possible to set a section to zero. We anticipate
that people might want to be able to come back from this by assigning a non-trivial
value to a section that was previously tuned to zero. This is why we keep such
trivial sections and do not delete them, say from `explicit_model_sections`
or `classes_of_model_sections`.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, defining_classes = Dict("b" => b))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(base_space(h)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice = Dict("b" => x2, "c0" => zero(parent(x1)))
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "b"  => x2
  "c0" => 0

julia> tuned_h = tune(h, my_choice)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(base_space(tuned_h)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice2 = Dict("b" => x2, "c0" => zero(parent(x1)))
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "b"  => x2
  "c0" => 0

julia> tuned_h2 = tune(tuned_h, my_choice2)
Hypersurface model over a concrete base

julia> is_zero(explicit_model_sections(tuned_h2)["c0"])
true

julia> x1, x2, x3 = gens(cox_ring(base_space(tuned_h2)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice3 = Dict("b" => x2, "c0" => x1^10)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "b"  => x2
  "c0" => x1^10

julia> tuned_h3 = tune(tuned_h2, my_choice3)
Hypersurface model over a concrete base

julia> is_zero(explicit_model_sections(tuned_h3)["c0"])
false
```
"""
function tune(h::HypersurfaceModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)
  # Consistency checks
  @req base_space(h) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return h
  secs_names = tunable_sections(h)
  tuned_secs_names = collect(keys(input_sections))
  @req all(in(secs_names), tuned_secs_names) "Provided section names are not among the tunable sections of the model"

  # 1. Tune model sections
  explicit_secs = deepcopy(explicit_model_sections(h))
  for x in tuned_secs_names
    section_parent = parent(input_sections[x])
    @req section_parent == parent(explicit_model_sections(h)[x]) "Parent mismatch between given and existing model section"
    if is_zero(input_sections[x]) == false
      @req degree(input_sections[x]) == divisor_class(classes_of_model_sections(h)[x]) "Degree mismatch between given and existing model section"
    end
    explicit_secs[x] = input_sections[x]
  end
  
  # 2. Compute the new hypersurface equation
  parametrized_hypersurface_equation = hypersurface_equation_parametrization(h)
  R = parent(parametrized_hypersurface_equation)
  vars = string.(symbols(R))
  S = cox_ring(ambient_space(h))
  images = [k in secs_names ? eval_poly(string(explicit_secs[k]), S) : k == "Kbar" ? eval_poly("0", S) : eval_poly(k, S) for k in vars]
  map = hom(R, S, images; check=false)
  new_hypersurface_equation = map(parametrized_hypersurface_equation)

  # 3. Build the new model
  model = HypersurfaceModel(explicit_secs, parametrized_hypersurface_equation, new_hypersurface_equation, base_space(h), ambient_space(h), fiber_ambient_space(h))
  set_attribute!(model, :partially_resolved, false)

  # 4. Copy the classes of model sections, but only of those sections that are used!
  new_classes_of_model_sections = Dict{String, ToricDivisorClass}()
  for key in keys(explicit_model_sections(model))
    m = divisor_class(classes_of_model_sections(h)[key]).coeff
    @req nrows(m) == 1 "Encountered inconsistency"
    new_classes_of_model_sections[key] = toric_divisor_class(base_space(model), m[1, :])
  end
  set_attribute!(model, :classes_of_model_sections => new_classes_of_model_sections)

  # 5. Return the model
  return model
end
