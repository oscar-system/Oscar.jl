#####################################################
# Setters
#####################################################

@doc raw"""
    set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)

Assigns a Weierstrass model to the given hypersurface model.

In the example below, we construct a hypersurface model and its corresponding 
Weierstrass model (see [BMT25](@cite BMT25) for background), and demonstrate how 
to associate them using this function.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3;

julia> h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p, completeness_check = false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(B2, weier_f, weier_g; completeness_check = false)
Weierstrass model over a concrete base

julia> set_weierstrass_model(h, w)
```
"""
function set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)
  set_attribute!(h, :weierstrass_model => w)
end

@doc raw"""
    set_global_tate_model(h::HypersurfaceModel, w::GlobalTateModel)

Assigns a global Tate model to the given hypersurface model.

In the example below, we construct a hypersurface model and its corresponding 
global Tate model (see [BMT25](@cite BMT25) for background), and demonstrate how 
to associate them using this function.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3;

julia> h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p, completeness_check = false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> a1 = 13 * x3^3;

julia> a2 = 7 * x1 * x2^5;

julia> a3 = x1^2 * x2^4 * x3^3;

julia> a4 = x1^3 * (x2 + x3)^9;

julia> a6 = zero(cox_ring(B2));

julia> t = global_tate_model(B2, [a1, a2, a3, a4, a6])
Global Tate model over a concrete base

julia> set_global_tate_model(h, t)
```
"""
function set_global_tate_model(h::HypersurfaceModel, t::GlobalTateModel)
  set_attribute!(h, :global_tate_model => t)
end



#####################################################
# 2: Tune a Hypersurface model
#####################################################

@doc raw"""
    tune(h::HypersurfaceModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)

Tune a hypersurface model by fixing a special choice for the model sections.
Note that it is in particular possible to set a section to zero. We anticipate
that people might want to be able to come back from this by assigning a non-trivial
value to a section that was previously tuned to zero. This is why we keep such
trivial sections and do not delete them, say from `explicit_model_sections`
or `classes_of_model_sections`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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
