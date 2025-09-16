@doc raw"""
    tune(w::WeierstrassModel, input_sections::Dict{String, <:Any})

Takes a Weierstrass model `w` and returns a new model in which the tunable
sections have been fixed to specific values provided by the user.

The `input_sections` argument is a dictionary mapping section names (as strings)
to values, typically given as elements of a multivariate polynomial ring. It is
also possible to set a section to zero.

Importantly, even if a section is tuned to zero, it is **not removed** from the
model’s metadata (such as `explicit_model_sections` or `classes_of_model_sections`).
This ensures the ability to later reintroduce non-trivial values for those sections,
preserving model flexibility and reversibility.

!!! note "Complete toric base"
    This function assumes that the toric base space is **complete**.
    Checking completeness may take a long time. To skip this check,
    pass the **optional keyword argument** `completeness_check=false`.

# Examples
```jldoctest
julia> using Random;

julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false, rng = Random.Xoshiro(1234))
Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> x1, x2, x3 = gens(coordinate_ring(base_space(w)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice = Dict("b" => x2, "c2" => x1^6)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "c2" => x1^6
  "b"  => x2

julia> tuned_w = tune(w, my_choice)
Weierstrass model over a concrete base

julia> explicit_model_sections(tuned_w)["c2"]
x1^6

julia> x1, x2, x3 = gens(coordinate_ring(base_space(tuned_w)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice2 = Dict("b" => x2, "c2" => zero(parent(x1)))
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "c2" => 0
  "b"  => x2

julia> tuned_w2 = tune(tuned_w, my_choice2)
Weierstrass model over a concrete base

julia> is_zero(explicit_model_sections(tuned_w2)["c2"])
true

julia> x1, x2, x3 = gens(coordinate_ring(base_space(tuned_w2)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice3 = Dict("b" => x2, "c2" => x1^6)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "c2" => x1^6
  "b"  => x2

julia> tuned_w3 = tune(tuned_w2, my_choice3)
Weierstrass model over a concrete base

julia> is_zero(explicit_model_sections(tuned_w3)["c2"])
false
```
"""
function tune(
  w::WeierstrassModel, input_sections::Dict{String,<:Any}; completeness_check::Bool=true
)
  # Consistency checks
  @req base_space(w) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return w
  secs_names = tunable_sections(w)
  tuned_secs_names = collect(keys(input_sections))
  @req all(in(secs_names), tuned_secs_names) "Provided section names are not among the tunable sections of the model"

  # 0. Prepare for computation by setting up some information
  explicit_secs = deepcopy(explicit_model_sections(w))
  def_secs_param = deepcopy(model_section_parametrization(w))
  weierstrass_sections = ["f", "g"]

  # 1. Tune model sections different from Weierstrass sections
  for x in setdiff(tuned_secs_names, weierstrass_sections)
    @req parent(input_sections[x]) == parent(explicit_model_sections(w)[x]) "Parent mismatch between given and existing model section"
    if is_zero(input_sections[x]) == false
      @req degree(input_sections[x]) == divisor_class(classes_of_model_sections(w)[x]) "Degree mismatch between given and existing model section"
    end
    explicit_secs[x] = input_sections[x]
  end

  # 2. Use model sections to reevaluate the Weierstrass sections via their known parametrization
  parametrization_keys = collect(keys(def_secs_param))
  if !isempty(parametrization_keys) && !isempty(secs_names)
    R = parent(def_secs_param[parametrization_keys[1]])
    S = parent(explicit_secs[secs_names[1]])
    vars = [string(k) for k in symbols(R)]
    images = [
      if k in secs_names
        explicit_secs[k]
      elseif k == "Kbar"
        eval_poly("0", S)
      else
        eval_poly(k, S)
      end for k in vars
    ]
    map = hom(R, S, images)
    for section in weierstrass_sections
      haskey(def_secs_param, section) &&
        (explicit_secs[section] = map(eval_poly(string(def_secs_param[section]), R)))
    end
  end

  # 3. Does the user want to set some Weierstrass sections? If so, overwrite existing choice with desired value.
  for sec in weierstrass_sections
    if haskey(input_sections, sec)
      @req parent(input_sections[sec]) == parent(explicit_model_sections(w)[sec]) "Parent mismatch between given and existing Weierstrass section"
      if is_zero(input_sections[sec]) == false
        @req degree(input_sections[sec]) == divisor_class(classes_of_model_sections(w)[sec]) "Degree mismatch between given and existing Weierstrass section"
      end
      explicit_secs[sec] = input_sections[sec]
      delete!(def_secs_param, sec)
    end
  end

  # 4. There could be unused model sections...
  if !isempty(parametrization_keys)
    polys = [
      eval_poly(string(def_secs_param[section]), R) for
      section in weierstrass_sections if haskey(def_secs_param, section)
    ]
    all_appearing_monomials = vcat([collect(monomials(p)) for p in polys]...)
    all_appearing_exponents = [collect(exponents(m))[1] for m in all_appearing_monomials]
    potentially_redundant_sections = gens(R)
    for k in 1:length(potentially_redundant_sections)
      string(potentially_redundant_sections[k]) in weierstrass_sections && continue
      is_used = any(
        all_appearing_exponents[l][k] != 0 for l in 1:length(all_appearing_exponents)
      )
      is_used || delete!(explicit_secs, string(potentially_redundant_sections[k]))
    end
  end

  # 5. After removing some sections, we must go over the parametrization again and adjust the ring in which the parametrization is given.
  if !isempty(def_secs_param)
    naive_vars = string.(gens(parent(first(values(def_secs_param)))))
    filtered_vars = filter(x -> haskey(explicit_secs, x), naive_vars)
    desired_ring, _ = polynomial_ring(QQ, filtered_vars; cached=false)
    for (key, value) in def_secs_param
      def_secs_param[key] = eval_poly(string(value), desired_ring)
    end
  end

  # 6. Build the new model
  resulting_model = weierstrass_model(
    base_space(w), explicit_secs, def_secs_param; completeness_check
  )

  # 7. Copy the classes of model sections
  new_classes_of_model_sections = Dict{String,ToricDivisorClass}()
  for (key, value) in classes_of_model_sections(w)
    m = divisor_class(value).coeff
    @req nrows(m) == 1 "Encountered inconsistency"
    new_classes_of_model_sections[key] = toric_divisor_class(
      base_space(resulting_model), m[1, :]
    )
  end
  set_attribute!(
    resulting_model, :classes_of_model_sections => new_classes_of_model_sections
  )

  # 8. Return the model
  return resulting_model
end

@doc raw"""
    tune(t::GlobalTateModel, input_sections::Dict{String, <:Any})

Tunes a global Tate model by specifying values for some of its defining sections.

This mechanism allows for specialized constructions—for example, setting certain sections to zero
to engineer enhanced singularities. Note that trivial (zero) sections are retained internally,
rather than deleted from attributes like `explicit_model_sections` or `classes_of_model_sections`.
This design choice enables users to later reintroduce nontrivial values for such sections
without loss of structure.

!!! note "Complete toric base"
    This function assumes that the toric base space is **complete**.
    Checking completeness may take a long time. To skip this check,
    pass the **optional keyword argument** `completeness_check=false`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> using Random;

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4 = gens(coordinate_ring(base_space(t)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice = Dict("a1" => x1^4, "w" => x2 - x3)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "w"  => x2 - x3
  "a1" => x1^4

julia> tuned_t = tune(t, my_choice)
Global Tate model over a concrete base

julia> tate_section_a1(tuned_t) == x1^4
true

julia> x1, x2, x3, x4 = gens(coordinate_ring(base_space(tuned_t)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice2 = Dict("a1" => zero(parent(x1)), "w" => x2 - x3)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "w"  => x2 - x3
  "a1" => 0

julia> tuned_t2 = tune(tuned_t, my_choice2)
Global Tate model over a concrete base

julia> is_zero(explicit_model_sections(tuned_t2)["a1"])
true

julia> x1, x2, x3, x4 = gens(coordinate_ring(base_space(tuned_t2)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice3 = Dict("a1" => x1^4, "w" => x2 - x3)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "w"  => x2 - x3
  "a1" => x1^4

julia> tuned_t3 = tune(tuned_t2, my_choice3)
Global Tate model over a concrete base

julia> is_zero(explicit_model_sections(tuned_t3)["a1"])
false
```
"""
function tune(
  t::GlobalTateModel, input_sections::Dict{String,<:Any}; completeness_check::Bool=true
)
  # Consistency checks
  @req base_space(t) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return t
  secs_names = tunable_sections(t)
  tuned_secs_names = collect(keys(input_sections))
  @req all(in(secs_names), tuned_secs_names) "Provided section names are not among the tunable sections of the model"

  # 0. Prepare for computation by setting up some information
  explicit_secs = deepcopy(explicit_model_sections(t))
  def_secs_param = deepcopy(model_section_parametrization(t))
  tate_sections = ["a1", "a2", "a3", "a4", "a6"]

  # 1. Tune model sections different from Tate sections
  for x in setdiff(tuned_secs_names, tate_sections)
    @req parent(input_sections[x]) == parent(explicit_model_sections(t)[x]) "Parent mismatch between given and existing model section"
    if is_zero(input_sections[x]) == false
      @req degree(input_sections[x]) == divisor_class(classes_of_model_sections(t)[x]) "Degree mismatch between given and existing model section"
    end
    explicit_secs[x] = input_sections[x]
  end

  # 2. Use model sections to reevaluate the Tate sections via their known parametrization
  parametrization_keys = collect(keys(def_secs_param))
  if !isempty(parametrization_keys) && !isempty(secs_names)
    R = parent(def_secs_param[parametrization_keys[1]])
    S = parent(explicit_secs[secs_names[1]])
    vars = [string(k) for k in symbols(R)]
    images = [
      if k in secs_names
        explicit_secs[k]
      elseif k == "Kbar"
        eval_poly("0", S)
      else
        eval_poly(k, S)
      end for k in vars
    ]
    map = hom(R, S, images)
    for section in tate_sections
      haskey(def_secs_param, section) &&
        (explicit_secs[section] = map(eval_poly(string(def_secs_param[section]), R)))
    end
  end

  # 3. Does the user want to set some Tate sections? If so, overwrite existing choice with desired value.
  for sec in tate_sections
    if haskey(input_sections, sec)
      @req parent(input_sections[sec]) == parent(explicit_model_sections(t)[sec]) "Parent mismatch between given and existing Tate section"
      if is_zero(input_sections[sec]) == false
        @req degree(input_sections[sec]) == divisor_class(classes_of_model_sections(t)[sec]) "Degree mismatch between given and existing Tate section"
      end
      explicit_secs[sec] = input_sections[sec]
      delete!(def_secs_param, sec)
    end
  end

  # 4. There could be unused model sections...
  if !isempty(parametrization_keys)
    polys = [
      eval_poly(string(def_secs_param[section]), R) for
      section in tate_sections if haskey(def_secs_param, section)
    ]
    all_appearing_monomials = vcat([collect(monomials(p)) for p in polys]...)
    all_appearing_exponents = [collect(exponents(m))[1] for m in all_appearing_monomials]
    potentially_redundant_sections = gens(R)
    for k in 1:length(potentially_redundant_sections)
      string(potentially_redundant_sections[k]) in tate_sections && continue
      is_used = any(
        all_appearing_exponents[l][k] != 0 for l in 1:length(all_appearing_exponents)
      )
      is_used || delete!(explicit_secs, string(potentially_redundant_sections[k]))
    end
  end

  # 5. After removing some sections, we must go over the parametrization again and adjust the ring in which the parametrization is given.
  if !isempty(def_secs_param)
    naive_vars = string.(gens(parent(first(values(def_secs_param)))))
    filtered_vars = filter(x -> haskey(explicit_secs, x), naive_vars)
    desired_ring, _ = polynomial_ring(QQ, filtered_vars; cached=false)
    for (key, value) in def_secs_param
      def_secs_param[key] = eval_poly(string(value), desired_ring)
    end
  end

  # 6. Build the new model
  resulting_model = global_tate_model(
    base_space(t), explicit_secs, def_secs_param; completeness_check
  )

  # 7. Copy the classes of model sections, but only of those sections that are used!
  new_classes_of_model_sections = Dict{String,ToricDivisorClass}()
  for key in keys(explicit_model_sections(resulting_model))
    m = divisor_class(classes_of_model_sections(t)[key]).coeff
    @req nrows(m) == 1 "Encountered inconsistency"
    new_classes_of_model_sections[key] = toric_divisor_class(
      base_space(resulting_model), m[1, :]
    )
  end
  set_attribute!(
    resulting_model, :classes_of_model_sections => new_classes_of_model_sections
  )

  # 8. Return the model
  return resulting_model
end

@doc raw"""
    tune(h::HypersurfaceModel, input_sections::Dict{String, <:Any})

Tune a hypersurface model by fixing a special choice for the model sections.
Note that it is in particular possible to set a section to zero. We anticipate
that people might want to be able to come back from this by assigning a non-trivial
value to a section that was previously tuned to zero. This is why we keep such
trivial sections and do not delete them, say from `explicit_model_sections`
or `classes_of_model_sections`.

!!! note "Complete toric base"
    This function assumes that the toric base space is **complete**.
    Checking completeness may take a long time. To skip this check,
    pass the **optional keyword argument** `completeness_check=false`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> using Random;

julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, defining_classes = Dict("b" => b), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(coordinate_ring(base_space(h)))
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

julia> x1, x2, x3 = gens(coordinate_ring(base_space(tuned_h)))
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

julia> x1, x2, x3 = gens(coordinate_ring(base_space(tuned_h2)))
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
function tune(
  h::HypersurfaceModel, input_sections::Dict{String,<:Any}; completeness_check::Bool=true
)
  # Consistency checks
  @req base_space(h) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return h
  secs_names = tunable_sections(h)
  tuned_secs_names = collect(keys(input_sections))
  @req all(in(secs_names), tuned_secs_names) "Provided section names are not among the tunable sections of the model"
  if completeness_check
    @req is_complete(base_space(h)) "Base space must be complete"
  end

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
  S = coordinate_ring(ambient_space(h))
  images = [
    if k in secs_names
      eval_poly(string(explicit_secs[k]), S)
    elseif k == "Kbar"
      eval_poly("0", S)
    else
      eval_poly(k, S)
    end for k in vars
  ]
  map = hom(R, S, images; check=false)
  new_hypersurface_equation = map(parametrized_hypersurface_equation)

  # 3. Build the new model
  model = HypersurfaceModel(
    explicit_secs,
    parametrized_hypersurface_equation,
    new_hypersurface_equation,
    base_space(h),
    ambient_space(h),
    fiber_ambient_space(h),
  )
  set_attribute!(model, :partially_resolved, false)

  # 4. Copy the classes of model sections, but only of those sections that are used!
  new_classes_of_model_sections = Dict{String,ToricDivisorClass}()
  for key in keys(explicit_model_sections(model))
    m = divisor_class(classes_of_model_sections(h)[key]).coeff
    @req nrows(m) == 1 "Encountered inconsistency"
    new_classes_of_model_sections[key] = toric_divisor_class(base_space(model), m[1, :])
  end
  set_attribute!(model, :classes_of_model_sections => new_classes_of_model_sections)

  # 5. Return the model
  return model
end

# FIXME: The below tune function is not consistent with our other "tune" functions (it does not correspond mathematically to a tuning)
# and produces an undesired form for explicit_model_sections. To be improved at a later date.

# @doc raw"""
#     tune(m::AbstractFTheoryModel, p::MPolyRingElem; completeness_check::Bool = true)

# Tune an F-theory model by replacing the hypersurface equation by a custom (polynomial)
# equation. The latter can be any type of polynomial: a Tate polynomial, a Weierstrass
# polynomial or a general polynomial. We do not conduct checks to tell which type the
# provided polynomial is. Consequently, this tuning will always return a hypersurface model.

# Note that there is less functionality for hypersurface models than for Weierstrass or Tate
# models. For instance, `singular_loci` can (currently) not be computed for hypersurface models.

# # Examples
# ```jldoctest
# julia> using Random;

# julia> B3 = projective_space(NormalToricVariety, 3)
# Normal toric variety

# julia> w = torusinvariant_prime_divisors(B3)[1]
# Torus-invariant, prime divisor on a normal toric variety

# julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false, rng = Random.Xoshiro(1234))
# Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

# julia> x1, x2, x3, x4, x, y, z = gens(parent(tate_polynomial(t)))
# 7-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
#  x1
#  x2
#  x3
#  x4
#  x
#  y
#  z

# julia> new_tate_polynomial = x^3 - y^2 - x * y * z * x4^4
# -x4^4*x*y*z + x^3 - y^2

# julia> tuned_t = tune(t, new_tate_polynomial)
# Hypersurface model over a concrete base

# julia> hypersurface_equation(tuned_t) == new_tate_polynomial
# true

# julia> base_space(tuned_t) == base_space(t)
# true
# ```
# """
# function tune(m::AbstractFTheoryModel, p::MPolyRingElem; completeness_check::Bool = true)
#   entry_test = (m isa GlobalTateModel) || (m isa WeierstrassModel) || (m isa HypersurfaceModel)
#   @req entry_test "Tuning currently supported only for Weierstrass, Tate and hypersurface models"
#   @req (base_space(m) isa NormalToricVariety) "Currently, tuning is only supported for models over concrete toric bases"
#   equation = hypersurface_equation(m)
#   @req parent(p) == parent(equation) "Parent mismatch between given and existing hypersurface polynomial"
#   @req degree(p) == degree(equation) "Degree mismatch between given and existing hypersurface polynomial"
#   p == equation && return m
#   explicit_model_sections = Dict{String, MPolyRingElem}()
#   gens_S = gens(parent(p))
#   for k in 1:length(gens_S)
#     explicit_model_sections[string(gens_S[k])] = gens_S[k]
#   end
#   tuned_model = HypersurfaceModel(explicit_model_sections, p, p, base_space(m), ambient_space(m), fiber_ambient_space(m))
#   set_attribute!(tuned_model, :partially_resolved, false)
#   return tuned_model
# end
