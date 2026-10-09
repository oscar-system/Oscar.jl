###############################################################################
#
#  Type-level data for reflections and reflecting hyperplanes
#
###############################################################################

"""
    ReflectionHyperplaneOrbitType

Describe one orbit of reflecting hyperplanes of a complex reflection group
type.

The descriptor records the component containing the orbit, its position among
the hyperplane orbits of that component, the order of the pointwise stabilizer
of a hyperplane, and the number of hyperplanes in the orbit. It also contains a
conventional representative word in the generator alphabet of the CHEVIE
reference model. Positive entries in the word denote generators and negative
entries denote their inverses.

The representative word and the orbit numbering are conventions attached to
the reference model; they are not intrinsic labels of the unmarked reflection
group type.
"""
struct ReflectionHyperplaneOrbitType
  group_type::ComplexReflectionGroupType
  component_index::Int
  local_orbit_index::Int
  pointwise_stabilizer_order::Int
  number_of_hyperplanes::ZZRingElem
  representative_word::Tuple{Vararg{Int}}
end

function Base.:(==)(
    O1::ReflectionHyperplaneOrbitType, O2::ReflectionHyperplaneOrbitType
  )
  return O1.group_type == O2.group_type &&
         O1.component_index == O2.component_index &&
         O1.local_orbit_index == O2.local_orbit_index &&
         O1.pointwise_stabilizer_order == O2.pointwise_stabilizer_order &&
         O1.number_of_hyperplanes == O2.number_of_hyperplanes &&
         O1.representative_word == O2.representative_word
end

function Base.isequal(
    O1::ReflectionHyperplaneOrbitType, O2::ReflectionHyperplaneOrbitType
  )
  return isequal(
    (
      O1.group_type,
      O1.component_index,
      O1.local_orbit_index,
      O1.pointwise_stabilizer_order,
      O1.number_of_hyperplanes,
      O1.representative_word,
    ),
    (
      O2.group_type,
      O2.component_index,
      O2.local_orbit_index,
      O2.pointwise_stabilizer_order,
      O2.number_of_hyperplanes,
      O2.representative_word,
    ),
  )
end

function Base.hash(O::ReflectionHyperplaneOrbitType, h::UInt)
  data = (
    O.group_type,
    O.component_index,
    O.local_orbit_index,
    O.pointwise_stabilizer_order,
    O.number_of_hyperplanes,
    O.representative_word,
  )
  return hash(data, hash(:ReflectionHyperplaneOrbitType, h))
end

"""
    ReflectionClassType

Describe a conjugacy class of reflections by a reflecting-hyperplane orbit and
an exponent.

If `O` is the underlying hyperplane-orbit descriptor and `s` is its
conventional distinguished reflection, then the descriptor with exponent `j`
represents the conjugacy class of `s^j`. The possible exponents are
`1:pointwise_stabilizer_order(O)-1`.
"""
struct ReflectionClassType
  orbit::ReflectionHyperplaneOrbitType
  exponent::Int

  function ReflectionClassType(
      orbit::ReflectionHyperplaneOrbitType, exponent::Int
    )
    @req 1 <= exponent < pointwise_stabilizer_order(orbit) "exponent must be between 1 and e_H - 1"
    return new(orbit, exponent)
  end
end

function Base.:(==)(C1::ReflectionClassType, C2::ReflectionClassType)
  return C1.orbit == C2.orbit && C1.exponent == C2.exponent
end

function Base.isequal(C1::ReflectionClassType, C2::ReflectionClassType)
  return isequal((C1.orbit, C1.exponent), (C2.orbit, C2.exponent))
end

function Base.hash(C::ReflectionClassType, h::UInt)
  return hash((C.orbit, C.exponent), hash(:ReflectionClassType, h))
end

"""
    complex_reflection_group_type(O::ReflectionHyperplaneOrbitType)
    complex_reflection_group_type(C::ReflectionClassType)

Return the complex reflection group type containing `O` or `C`.
"""
complex_reflection_group_type(O::ReflectionHyperplaneOrbitType) = O.group_type
complex_reflection_group_type(C::ReflectionClassType) =
  complex_reflection_group_type(hyperplane_orbit(C))

"""
    component_index(O::ReflectionHyperplaneOrbitType)
    component_index(C::ReflectionClassType)

Return the index of the irreducible component containing `O` or `C`.
"""
component_index(O::ReflectionHyperplaneOrbitType) = O.component_index
component_index(C::ReflectionClassType) = component_index(hyperplane_orbit(C))

"""
    local_orbit_index(O::ReflectionHyperplaneOrbitType)
    local_orbit_index(C::ReflectionClassType)

Return the position of the hyperplane orbit within its irreducible component.
"""
local_orbit_index(O::ReflectionHyperplaneOrbitType) = O.local_orbit_index
local_orbit_index(C::ReflectionClassType) = local_orbit_index(hyperplane_orbit(C))

"""
    pointwise_stabilizer_order(O::ReflectionHyperplaneOrbitType)
    pointwise_stabilizer_order(C::ReflectionClassType)

Return the order `e_H` of the pointwise stabilizer of a hyperplane in the
underlying orbit.
"""
pointwise_stabilizer_order(O::ReflectionHyperplaneOrbitType) =
  O.pointwise_stabilizer_order
pointwise_stabilizer_order(C::ReflectionClassType) =
  pointwise_stabilizer_order(hyperplane_orbit(C))

"""
    number_of_hyperplanes(O::ReflectionHyperplaneOrbitType)
    number_of_hyperplanes(C::ReflectionClassType)

Return the number of hyperplanes in the underlying orbit. This is also the
number of reflections in the represented conjugacy class.
"""
number_of_hyperplanes(O::ReflectionHyperplaneOrbitType) = O.number_of_hyperplanes
number_of_hyperplanes(C::ReflectionClassType) =
  number_of_hyperplanes(hyperplane_orbit(C))

"""
    orbit_size(O::ReflectionHyperplaneOrbitType)
    orbit_size(C::ReflectionClassType)

Return the number of hyperplanes in the underlying hyperplane orbit.
"""
orbit_size(O::ReflectionHyperplaneOrbitType) = number_of_hyperplanes(O)
orbit_size(C::ReflectionClassType) = number_of_hyperplanes(C)

"""
    hyperplane_orbit(C::ReflectionClassType)

Return the reflecting-hyperplane orbit underlying `C`.
"""
hyperplane_orbit(C::ReflectionClassType) = C.orbit

"""
    reflection_exponent(C::ReflectionClassType)

Return the exponent of the conventional distinguished reflection represented
by `C`.
"""
reflection_exponent(C::ReflectionClassType) = C.exponent

"""
    representative_word(O::ReflectionHyperplaneOrbitType)
    representative_word(C::ReflectionClassType)

Return a conventional representative word in the generator alphabet of the
CHEVIE reference model.
"""
representative_word(O::ReflectionHyperplaneOrbitType) = O.representative_word
function representative_word(C::ReflectionClassType)
  word = representative_word(hyperplane_orbit(C))
  return Tuple(x for _ in 1:reflection_exponent(C) for x in word)
end

function Base.show(io::IO, O::ReflectionHyperplaneOrbitType)
  print(
    io,
    "Reflection hyperplane orbit ",
    local_orbit_index(O),
    " in component ",
    component_index(O),
    " (e_H = ",
    pointwise_stabilizer_order(O),
    ", ",
    number_of_hyperplanes(O),
    " hyperplanes)",
  )
end

function Base.show(io::IO, C::ReflectionClassType)
  print(
    io,
    "Reflection class in component ",
    component_index(C),
    ", hyperplane orbit ",
    local_orbit_index(C),
    ", exponent ",
    reflection_exponent(C),
  )
end

###############################################################################
#
#  Irreducible exceptional types
#
###############################################################################

# Each entry is `(reference generator, e_H, number of hyperplanes)`. The outer
# tuple is indexed by the Shephard--Todd number minus three. Within each entry,
# the orbits are ordered by the index of their CHEVIE reference generator. The
# totals agree with Lehrer--Taylor (2009), Table D.3, and the orbit partitions
# below were validated against every available CHEVIE and Magma model (and
# every available Lehrer--Taylor model).
const _exceptional_reflection_hyperplane_orbit_table = (
  ((1, 3, 4),),                                      # G4
  ((1, 3, 4), (2, 3, 4)),                            # G5
  ((1, 2, 6), (2, 3, 4)),                            # G6
  ((1, 2, 6), (2, 3, 4), (3, 3, 4)),                # G7
  ((1, 4, 6),),                                      # G8
  ((1, 2, 12), (2, 4, 6)),                           # G9
  ((1, 3, 8), (2, 4, 6)),                            # G10
  ((1, 2, 12), (2, 3, 8), (3, 4, 6)),               # G11
  ((1, 2, 12),),                                     # G12
  ((1, 2, 6), (2, 2, 12)),                           # G13
  ((1, 2, 12), (2, 3, 8)),                           # G14
  ((1, 2, 12), (2, 3, 8), (3, 2, 6)),               # G15
  ((1, 5, 12),),                                     # G16
  ((1, 2, 30), (2, 5, 12)),                          # G17
  ((1, 3, 20), (2, 5, 12)),                          # G18
  ((1, 2, 30), (2, 3, 20), (3, 5, 12)),             # G19
  ((1, 3, 20),),                                     # G20
  ((1, 2, 30), (2, 3, 20)),                          # G21
  ((1, 2, 30),),                                     # G22
  ((1, 2, 15),),                                     # G23
  ((1, 2, 21),),                                     # G24
  ((1, 3, 12),),                                     # G25
  ((1, 2, 9), (2, 3, 12)),                           # G26
  ((1, 2, 45),),                                     # G27
  ((1, 2, 12), (3, 2, 12)),                          # G28
  ((1, 2, 40),),                                     # G29
  ((1, 2, 60),),                                     # G30
  ((1, 2, 60),),                                     # G31
  ((1, 3, 40),),                                     # G32
  ((1, 2, 45),),                                     # G33
  ((1, 2, 126),),                                    # G34
  ((1, 2, 36),),                                     # G35
  ((1, 2, 63),),                                     # G36
  ((1, 2, 120),),                                    # G37
)

function _exceptional_reflection_hyperplane_orbits(
    G::ComplexReflectionGroupType, component::Int, st_number::Int
  )
  data = _exceptional_reflection_hyperplane_orbit_table[st_number - 3]
  @assert issorted(first.(data))

  result = ReflectionHyperplaneOrbitType[]
  for (local_index, (generator, stabilizer_order, size)) in enumerate(data)
    push!(
      result,
      ReflectionHyperplaneOrbitType(
        G, component, local_index, stabilizer_order, ZZ(size), (generator,)
      ),
    )
  end
  return result
end

###############################################################################
#
#  The infinite series
#
###############################################################################

function _imprimitive_reflection_hyperplane_orbits(
    G::ComplexReflectionGroupType, component::Int, type::Tuple{Int, Int, Int}
  )
  m, p, n = type
  @assert is_divisible_by(m, p)

  result = ReflectionHyperplaneOrbitType[]

  # A normalized rank-one type is G(m, 1, 1). Its unique hyperplane is the
  # origin, and its pointwise stabilizer is the whole cyclic group.
  if n == 1
    @assert p == 1
    push!(result, ReflectionHyperplaneOrbitType(G, component, 1, m, ZZ(1), (1,)))
    return result
  end

  local_index = 0

  # The coordinate hyperplanes form one orbit when they occur. In the CHEVIE
  # reference generators, the first generator represents this orbit.
  if m > p
    local_index += 1
    push!(
      result,
      ReflectionHyperplaneOrbitType(
        G, component, local_index, div(m, p), ZZ(n), (1,)
      ),
    )
  end

  number_of_transposition_hyperplanes = div(ZZ(m) * n * (n - 1), 2)

  # If a coordinate-reflection generator occurs, the first transposition-type
  # generator has index two; otherwise it has index one. For G(1, 1, n), this
  # is the first simple reflection of the essential symmetric-group model.
  transposition_generator = m > p ? 2 : 1

  if n == 2 && is_even(p)
    # In this case the transposition hyperplanes split into two equally sized
    # orbits, represented by two consecutive CHEVIE reference generators.
    @assert is_even(number_of_transposition_hyperplanes)
    size = div(number_of_transposition_hyperplanes, 2)
    for generator in transposition_generator:(transposition_generator + 1)
      local_index += 1
      push!(
        result,
        ReflectionHyperplaneOrbitType(
          G, component, local_index, 2, size, (generator,)
        ),
      )
    end
  else
    local_index += 1
    push!(
      result,
      ReflectionHyperplaneOrbitType(
        G,
        component,
        local_index,
        2,
        number_of_transposition_hyperplanes,
        (transposition_generator,),
      ),
    )
  end

  return result
end

###############################################################################
#
#  Public catalogue functions
#
###############################################################################

"""
    reflection_hyperplane_orbits(G::ComplexReflectionGroupType)

Return type-level descriptors for the orbits of reflecting hyperplanes of `G`.

For a direct product, the result is the concatenation of the component
catalogues. Each descriptor retains its component index and a representative
word in that component's CHEVIE reference generator alphabet.
"""
function reflection_hyperplane_orbits(G::ComplexReflectionGroupType)
  result = ReflectionHyperplaneOrbitType[]

  for (component, type) in enumerate(G.type)
    if type isa Int
      append!(result, _exceptional_reflection_hyperplane_orbits(G, component, type))
    else
      append!(result, _imprimitive_reflection_hyperplane_orbits(G, component, type))
    end
  end

  @assert sum(number_of_hyperplanes, result; init=ZZ(0)) == number_of_hyperplanes(G)
  @assert sum(
    O -> number_of_hyperplanes(O) * (pointwise_stabilizer_order(O) - 1),
    result;
    init=ZZ(0),
  ) == number_of_reflections(G)

  return result
end

"""
    reflection_classes(O::ReflectionHyperplaneOrbitType)
    reflection_classes(G::ComplexReflectionGroupType)

Return type-level descriptors for reflection conjugacy classes.

The classes over a hyperplane orbit `O` are derived as the pairs `(O, j)` for
`1 <= j < pointwise_stabilizer_order(O)`. For a group type, the classes are
ordered first by component, then by the local hyperplane-orbit index, and then
by the exponent.
"""
function reflection_classes(O::ReflectionHyperplaneOrbitType)
  return [ReflectionClassType(O, j) for j in 1:(pointwise_stabilizer_order(O) - 1)]
end

function reflection_classes(G::ComplexReflectionGroupType)
  result = ReflectionClassType[]
  for O in reflection_hyperplane_orbits(G)
    append!(result, reflection_classes(O))
  end
  return result
end
