###############################################################################
#
#  Reflection data in a concrete matrix realization
#
###############################################################################

"""
    ComplexReflectionHyperplaneOrbit

An orbit of reflecting hyperplanes in a marked matrix realization of a complex
reflection group.

Besides the type-level descriptor, this object stores the concrete hyperplanes
and one distinguished generator of their cyclic pointwise stabilizer for every
hyperplane. The latter form a conjugacy orbit. Their powers give all
reflections over the hyperplane orbit.
"""
struct ComplexReflectionHyperplaneOrbit{G<:MatGroup, H<:GSet, R<:GSet}
  group::G
  orbit_type::ReflectionHyperplaneOrbitType
  hyperplanes::H
  distinguished_reflections::R
end

"""
    ComplexReflectionClass

A conjugacy class of reflections in a marked matrix realization.

The class is labelled by a [`ReflectionClassType`](@ref). Its elements are
computed as the conjugacy orbit of the corresponding power of a distinguished
reflection.
"""
struct ComplexReflectionClass{O<:ComplexReflectionHyperplaneOrbit, R<:GSet}
  class_type::ReflectionClassType
  orbit::O
  elements::R
end

function Base.show(io::IO, O::ComplexReflectionHyperplaneOrbit)
  print(
    io,
    "Concrete reflection hyperplane orbit ",
    local_orbit_index(O),
    " in component ",
    component_index(O),
    " with ",
    number_of_hyperplanes(O),
    " hyperplanes",
  )
end

function Base.show(io::IO, C::ComplexReflectionClass)
  print(
    io,
    "Concrete reflection class in component ",
    component_index(C),
    ", hyperplane orbit ",
    local_orbit_index(C),
    ", exponent ",
    reflection_exponent(C),
  )
end

###############################################################################
#
#  Getters
#
###############################################################################

parent(O::ComplexReflectionHyperplaneOrbit) = O.group
parent(C::ComplexReflectionClass) = parent(hyperplane_orbit(C))

"""
    reflection_hyperplane_orbit_type(O::ComplexReflectionHyperplaneOrbit)

Return the type-level reference descriptor underlying the concrete hyperplane
orbit `O`.
"""
reflection_hyperplane_orbit_type(O::ComplexReflectionHyperplaneOrbit) = O.orbit_type

"""
    reflection_class_type(C::ComplexReflectionClass)

Return the type-level reference descriptor underlying the concrete reflection
class `C`.
"""
reflection_class_type(C::ComplexReflectionClass) = C.class_type

complex_reflection_group_type(O::ComplexReflectionHyperplaneOrbit) =
  complex_reflection_group_type(reflection_hyperplane_orbit_type(O))
complex_reflection_group_type(C::ComplexReflectionClass) =
  complex_reflection_group_type(reflection_class_type(C))

component_index(O::ComplexReflectionHyperplaneOrbit) =
  component_index(reflection_hyperplane_orbit_type(O))
component_index(C::ComplexReflectionClass) = component_index(reflection_class_type(C))

local_orbit_index(O::ComplexReflectionHyperplaneOrbit) =
  local_orbit_index(reflection_hyperplane_orbit_type(O))
local_orbit_index(C::ComplexReflectionClass) = local_orbit_index(reflection_class_type(C))

pointwise_stabilizer_order(O::ComplexReflectionHyperplaneOrbit) =
  pointwise_stabilizer_order(reflection_hyperplane_orbit_type(O))
pointwise_stabilizer_order(C::ComplexReflectionClass) =
  pointwise_stabilizer_order(reflection_class_type(C))

number_of_hyperplanes(O::ComplexReflectionHyperplaneOrbit) = length(O.hyperplanes)
number_of_hyperplanes(C::ComplexReflectionClass) =
  number_of_hyperplanes(hyperplane_orbit(C))

orbit_size(O::ComplexReflectionHyperplaneOrbit) = number_of_hyperplanes(O)
orbit_size(C::ComplexReflectionClass) = length(C.elements)

representative_word(O::ComplexReflectionHyperplaneOrbit) =
  representative_word(reflection_hyperplane_orbit_type(O))
representative_word(C::ComplexReflectionClass) = representative_word(reflection_class_type(C))

reflection_exponent(C::ComplexReflectionClass) =
  reflection_exponent(reflection_class_type(C))

hyperplane_orbit(C::ComplexReflectionClass) = C.orbit

"""
    reflection_hyperplanes(O::ComplexReflectionHyperplaneOrbit)

Return the concrete reflecting hyperplanes in `O` as a `GSet`.
"""
reflection_hyperplanes(O::ComplexReflectionHyperplaneOrbit) = O.hyperplanes

"""
    distinguished_reflections(O::ComplexReflectionHyperplaneOrbit)

Return the conjugacy orbit of the marked generators of the cyclic pointwise
stabilizers of the hyperplanes in `O`.
"""
distinguished_reflections(O::ComplexReflectionHyperplaneOrbit) =
  O.distinguished_reflections
elements(C::ComplexReflectionClass) = C.elements

representative(O::ComplexReflectionHyperplaneOrbit) =
  representative(distinguished_reflections(O))
representative(C::ComplexReflectionClass) = representative(elements(C))

Base.length(O::ComplexReflectionHyperplaneOrbit) = length(distinguished_reflections(O))
Base.length(C::ComplexReflectionClass) = length(elements(C))

###############################################################################
#
#  Construction
#
###############################################################################

_conjugate_reflection(s::MatGroupElem, g::MatGroupElem) = inv(g) * s * g

function _fixed_hyperplane(s::MatGroupElem)
  is_reflection, data = is_complex_reflection_with_data(matrix(s))
  @assert is_reflection
  return hyperplane(data)
end

function _distinguished_reflection_pairs(O::ComplexReflectionHyperplaneOrbit)
  return [(_fixed_hyperplane(s), s) for s in distinguished_reflections(O)]
end

"""
    distinguished_reflection(O::ComplexReflectionHyperplaneOrbit, H)

Return the distinguished generator of the cyclic pointwise stabilizer of the
reflecting hyperplane `H` in `O`.
"""
function distinguished_reflection(O::ComplexReflectionHyperplaneOrbit, H)
  @req H in reflection_hyperplanes(O) "The hyperplane does not belong to this orbit"
  pairs = _distinguished_reflection_pairs(O)
  position = findfirst(pair -> first(pair) == H, pairs)
  @assert position !== nothing
  return last(pairs[position])
end

function _reflection_hyperplane_orbit_candidates(W::MatGroup, seeds=gens(W))
  candidates = []

  for s in seeds
    is_reflection, data = is_complex_reflection_with_data(s)
    @req is_reflection "The marked generators must be complex reflections"

    H = hyperplane(data)
    if any(candidate -> H in candidate.hyperplanes, candidates)
      continue
    end

    hyperplane_orbit = gset(W, [H])
    push!(
      candidates,
      (
        order=Int(order(data)),
        hyperplanes=hyperplane_orbit,
        seed=s,
      ),
    )
  end

  return candidates
end

function _imprimitive_reference_generator_indices(
    type::Tuple{Int, Int, Int}, model::Symbol
  )
  m, p, n = type
  n == 1 && return [1]

  split_transposition_orbit = n == 2 && is_even(p)
  if model == :CHEVIE
    if m == p
      return split_transposition_orbit ? [1, 2] : [1]
    elseif p == 1
      return [1, 2]
    else
      return split_transposition_orbit ? [1, 2, 3] : [1, 2]
    end
  elseif model == :LT
    if m == p
      return split_transposition_orbit ? [1, 2] : [1]
    elseif p == 1
      return [1, 2]
    else
      return split_transposition_orbit ? [2, 1, 3] : [2, 1]
    end
  elseif model == :Magma
    if m == p
      return split_transposition_orbit ? [n, 1] : [1]
    elseif p == 1
      return [n, 1]
    else
      return split_transposition_orbit ? [n + 1, n, 1] : [n + 1, n]
    end
  end

  return nothing
end

# The correspondence between reference hyperplane orbits and the generators
# of every built-in exceptional model is marking data. It is recorded in full
# rather than inferred from `(e_H, orbit size)`, since those invariants do not
# distinguish all orbits. The tuples are indexed by the Shephard--Todd number
# minus three and then by the local reference-orbit index.
const _LT_EXCEPTIONAL_ORBIT_GENERATOR_INDICES = (
  (1,),       # G4
  (1, 2),     # G5
  (1, 2),     # G6
  (1, 2, 3),  # G7
  (1,),       # G8
  (1, 2),     # G9
  (1, 2),     # G10
  (2, 1, 3),  # G11
  (1,),       # G12
  (1, 2),     # G13
  (2, 1),     # G14
  (3, 2, 1),  # G15
  (1,),       # G16
  (2, 1),     # G17
  (1, 2),     # G18
  (1, 2, 3),  # G19
  (1,),       # G20
  (1, 2),     # G21
  (1,),       # G22
  (1,),       # G23
)

const _MAGMA_EXCEPTIONAL_ORBIT_GENERATOR_INDICES = (
  (1,),       # G4
  (1, 2),     # G5
  (1, 2),     # G6
  (1, 2, 3),  # G7
  (1,),       # G8
  (1, 2),     # G9
  (2, 1),     # G10
  (2, 3, 1),  # G11
  (1,),       # G12
  (1, 2),     # G13
  (2, 1),     # G14
  (3, 2, 1),  # G15
  (1,),       # G16
  (1, 2),     # G17
  (1, 2),     # G18
  (3, 2, 1),  # G19
  (1,),       # G20
  (1, 2),     # G21
  (1,),       # G22
  (1,),       # G23
  (1,),       # G24
  (1,),       # G25
  (3, 1),     # G26
  (1,),       # G27
  (1, 3),     # G28
  (1,),       # G29
  (1,),       # G30
  (1,),       # G31
  (1,),       # G32
  (1,),       # G33
  (1,),       # G34
  (1,),       # G35
  (1,),       # G36
  (1,),       # G37
)

function _exceptional_model_orbit_generator_indices(type::Int, model::Symbol)
  if model == :LT && 4 <= type <= 23
    return _LT_EXCEPTIONAL_ORBIT_GENERATOR_INDICES[type - 3]
  elseif model == :Magma && 4 <= type <= 37
    return _MAGMA_EXCEPTIONAL_ORBIT_GENERATOR_INDICES[type - 3]
  end
  return nothing
end

# A model generator can represent a proper power of the CHEVIE reference
# reflection, even after its hyperplane orbit has been identified. These are
# the powers that turn the selected model generator into the distinguished
# reflection of the reference marking. They are determined by comparing the
# unique nontrivial eigenvalues under the presentation embeddings used by the
# constructors: the named positive real roots and standard complex roots in
# the Lehrer--Taylor towers, and the standard cyclotomic generators in the
# CHEVIE and Magma models. Unlisted built-in cases use power one.
const _exceptional_model_orbit_reference_powers = Dict(
  (:LT, 4, 1) => 2,
  (:LT, 5, 1) => 2,
  (:LT, 5, 2) => 2,
  (:LT, 6, 2) => 2,
  (:LT, 7, 2) => 2,
  (:LT, 7, 3) => 2,
  (:LT, 10, 1) => 2,
  (:LT, 11, 2) => 2,
  (:LT, 14, 2) => 2,
  (:LT, 15, 2) => 2,
  (:LT, 16, 1) => 4,
  (:LT, 17, 2) => 4,
  (:LT, 18, 2) => 4,
  (:LT, 19, 2) => 2,
  (:LT, 19, 3) => 4,
  (:LT, 20, 1) => 2,
  (:LT, 21, 2) => 2,
  (:Magma, 9, 2) => 3,
  (:Magma, 10, 2) => 3,
  (:Magma, 11, 3) => 3,
)

function _marked_generator_index(
    W::MatGroup, orbit_type::ReflectionHyperplaneOrbitType
  )
  group_type = complex_reflection_group_type(W)
  component = component_index(orbit_type)
  type = group_type.type[component]
  models = complex_reflection_group_model(W)
  models isa AbstractVector || return nothing

  component_groups = components(W)
  component_groups === nothing && return nothing
  offset = sum(i -> length(gens(component_groups[i])), 1:(component - 1); init=0)

  if type isa Tuple
    indices = _imprimitive_reference_generator_indices(type, models[component])
    indices === nothing && return nothing
    return offset + indices[local_orbit_index(orbit_type)]
  elseif models[component] == :CHEVIE
    # Exceptional type-level words are defined in the CHEVIE marking.
    return offset + first(representative_word(orbit_type))
  else
    indices = _exceptional_model_orbit_generator_indices(
      type, models[component]
    )
    if indices !== nothing
      return offset + indices[local_orbit_index(orbit_type)]
    end
  end

  return nothing
end

function _model_reference_power(
    W::MatGroup, orbit_type::ReflectionHyperplaneOrbitType
  )
  group_type = complex_reflection_group_type(W)
  component = component_index(orbit_type)
  type = group_type.type[component]
  models = complex_reflection_group_model(W)
  models isa AbstractVector || return nothing
  model = models[component]
  model in (:CHEVIE, :LT, :Magma) || return nothing
  type isa Tuple && return 1
  return get(
    _exceptional_model_orbit_reference_powers,
    (model, type, local_orbit_index(orbit_type)),
    1,
  )
end

function _has_known_reflection_hyperplane_orbit_marking(W::MatGroup)
  return has_attribute(W, :reflection_hyperplane_orbit_marking_is_known) &&
         get_attribute(W, :reflection_hyperplane_orbit_marking_is_known)
end

"""
    set_reflection_hyperplane_orbit_marking!(W::MatGroup, representatives;
                                             class_exponents=nothing)

Mark a custom matrix realization by giving the image of the reference
distinguished reflection for each type-level reflecting-hyperplane orbit, in
the order returned by
`reflection_hyperplane_orbits(complex_reflection_group_type(W))`.

Each reflection must generate the cyclic pointwise stabilizer of its fixed
hyperplane. By default, each supplied reflection is understood to be the image
of the reference distinguished reflection itself. If a supplied reflection is
the `u`-th power of that reference, pass `u` in `class_exponents`; the values
must be invertible modulo the corresponding pointwise-stabilizer orders. The
stored representatives are normalized back to reference exponent one.

Built-in models already provide a marking compatible with the CHEVIE reference
words. The marking must be set before the concrete hyperplane orbits or
reflection classes are computed.
"""
function set_reflection_hyperplane_orbit_marking!(
    W::MatGroup,
    representatives::AbstractVector;
    class_exponents=nothing,
  )
  group_type = complex_reflection_group_type(W)
  @req group_type !== nothing "The matrix group has no assigned complex reflection group type"
  @req !has_attribute(
    W, :reflection_hyperplane_orbits
  ) "The concrete hyperplane orbits have already been computed"

  orbit_types = reflection_hyperplane_orbits(group_type)
  @req length(representatives) == length(
    orbit_types
  ) "Give one representative for each type-level hyperplane orbit"
  if class_exponents === nothing
    class_exponents = ones(Int, length(orbit_types))
  else
    @req length(class_exponents) == length(
      orbit_types
    ) "Give one class exponent for each marking representative"
  end

  marking = elem_type(W)[]
  for (s, u, orbit_type) in zip(representatives, class_exponents, orbit_types)
    @req s isa MatGroupElem && parent(
      s
    ) === W "Every representative must be an element of the matrix group"
    is_reflection, data = is_complex_reflection_with_data(s)
    @req is_reflection "Every marking representative must be a complex reflection"
    @req order(data) == pointwise_stabilizer_order(
      orbit_type
    ) "A marking representative does not generate the full pointwise stabilizer"
    @req u isa IntegerUnion "Every class exponent must be an integer"
    e = pointwise_stabilizer_order(orbit_type)
    normalized_exponent = Int(mod(u, e))
    @req gcd(normalized_exponent, e) == 1 "Every class exponent must be invertible modulo the pointwise-stabilizer order"
    push!(marking, s^invmod(normalized_exponent, e))
  end

  set_attribute!(W, :reflection_hyperplane_orbit_marking, marking)
  set_attribute!(W, :reflection_hyperplane_orbit_marking_is_known, true)
  return W
end

"""
    reflection_hyperplane_orbits(W::MatGroup)

Return the reflecting-hyperplane orbits of the marked complex reflection group
`W`.

This method uses the known Shephard--Todd type of `W`; it does not try to
recognize an arbitrary matrix group. Orbit labels come from the reference
catalogue, while hyperplanes and distinguished reflections are computed in the
given matrix realization. The computation uses matrix and subspace orbits in
Julia and does not require GAP conjugacy classes, so it also works for the
relative number-field towers used by the Lehrer--Taylor models.
"""
function reflection_hyperplane_orbits(W::MatGroup)
  if has_attribute(W, :reflection_hyperplane_orbits)
    return get_attribute(W, :reflection_hyperplane_orbits)
  end

  group_type = complex_reflection_group_type(W)
  @req group_type !== nothing "The matrix group has no assigned complex reflection group type"

  orbit_types = reflection_hyperplane_orbits(group_type)
  supplied_marking = has_attribute(W, :reflection_hyperplane_orbit_marking) ?
                     get_attribute(W, :reflection_hyperplane_orbit_marking) : nothing
  candidates = _reflection_hyperplane_orbit_candidates(
    W, supplied_marking === nothing ? gens(W) : supplied_marking
  )
  unused = trues(length(candidates))
  result = ComplexReflectionHyperplaneOrbit[]

  for (orbit_index, orbit_type) in enumerate(orbit_types)
    model_reflection = if supplied_marking !== nothing
      supplied_marking[orbit_index]
    else
      marked_generator_index = _marked_generator_index(W, orbit_type)
      marked_generator_index === nothing ? nothing : gens(W)[marked_generator_index]
    end

    position = if model_reflection === nothing
      positions = findall(eachindex(candidates)) do i
        unused[i] &&
          candidates[i].order == pointwise_stabilizer_order(orbit_type) &&
          length(candidates[i].hyperplanes) == number_of_hyperplanes(orbit_type)
      end
      @req length(positions) <= 1 "The type data do not distinguish several " *
                                  "hyperplane orbits; supply an explicit marking"
      isempty(positions) ? nothing : first(positions)
    else
      marked_hyperplane = _fixed_hyperplane(model_reflection)
      findfirst(eachindex(candidates)) do i
        unused[i] && marked_hyperplane in candidates[i].hyperplanes
      end
    end

    @req position !== nothing "The model generators do not realize the " *
                              "reference hyperplane-orbit marking"
    unused[position] = false
    candidate = candidates[position]
    @assert candidate.order == pointwise_stabilizer_order(orbit_type)
    @assert length(candidate.hyperplanes) == number_of_hyperplanes(orbit_type)

    if model_reflection === nothing
      model_reflection = candidate.seed
    end
    reference_reflection = if supplied_marking !== nothing
      model_reflection
    else
      reference_power = _model_reference_power(W, orbit_type)
      if reference_power === nothing
        @req pointwise_stabilizer_order(orbit_type) == 2 "The type and orbit " *
          "size do not determine the reflection-class orientation; supply an " *
          "explicit marking"
        model_reflection
      else
        model_reflection^reference_power
      end
    end
    reflection_orbit = gset(W, _conjugate_reflection, [reference_reflection])

    # The normalizer of a reflecting hyperplane centralizes its pointwise
    # stabilizer (Lehrer--Taylor (2009), Proposition 1.20). Thus conjugating a
    # distinguished reflection gives exactly one distinguished reflection for
    # each hyperplane in the orbit.
    @assert length(candidate.hyperplanes) == length(reflection_orbit)
    push!(
      result,
      ComplexReflectionHyperplaneOrbit(
        W, orbit_type, candidate.hyperplanes, reflection_orbit
      ),
    )
  end

  @assert sum(number_of_hyperplanes, result; init=0) == number_of_hyperplanes(group_type)
  set_attribute!(W, :reflection_hyperplane_orbits, result)
  set_attribute!(W, :reflection_hyperplane_orbit_marking_is_known, true)
  return result
end

"""
    reflection_classes(O::ComplexReflectionHyperplaneOrbit)
    reflection_classes(W::MatGroup)

Return the concrete conjugacy classes of reflections over `O`, or all concrete
reflection classes of the marked complex reflection group `W`.

For a hyperplane orbit with pointwise stabilizer of order `e_H`, the classes
are the conjugacy orbits of `s^j` for `1 <= j < e_H`, where `s` is its marked
distinguished reflection.
"""
function reflection_classes(O::ComplexReflectionHyperplaneOrbit)
  W = parent(O)
  s = representative(O)
  result = ComplexReflectionClass[]

  for class_type in reflection_classes(reflection_hyperplane_orbit_type(O))
    class_elements = gset(
      W,
      _conjugate_reflection,
      [s^reflection_exponent(class_type)],
    )
    @assert length(class_elements) == number_of_hyperplanes(O)
    push!(result, ComplexReflectionClass(class_type, O, class_elements))
  end

  return result
end

function reflection_classes(W::MatGroup)
  if has_attribute(W, :reflection_classes)
    return get_attribute(W, :reflection_classes)
  end

  result = ComplexReflectionClass[]
  for O in reflection_hyperplane_orbits(W)
    append!(result, reflection_classes(O))
  end

  @assert length(result) == number_of_reflection_classes(complex_reflection_group_type(W))
  set_attribute!(W, :reflection_classes, result)
  return result
end

"""
    reflection_library(O::ComplexReflectionHyperplaneOrbit)
    reflection_library(W::MatGroup)

Return the reflections arranged first by hyperplane orbit, then by hyperplane,
and finally by nonidentity element of the cyclic pointwise stabilizer.

For one orbit `O`, the result is a vector indexed by its hyperplanes. For `W`,
the result has one additional outer level indexed by hyperplane orbits. No
ordering of individual hyperplanes is part of the type-level marking.
"""
function reflection_library(O::ComplexReflectionHyperplaneOrbit)
  e = pointwise_stabilizer_order(O)
  pairs = _distinguished_reflection_pairs(O)
  return [
    begin
      position = findfirst(pair -> first(pair) == H, pairs)
      @assert position !== nothing
      s = last(pairs[position])
      [s^j for j in 1:(e - 1)]
    end for H in reflection_hyperplanes(O)
  ]
end

reflection_library(W::MatGroup) =
  [reflection_library(O) for O in reflection_hyperplane_orbits(W)]
