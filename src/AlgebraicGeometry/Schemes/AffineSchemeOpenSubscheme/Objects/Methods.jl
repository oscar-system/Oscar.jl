
########################################################################
# Implementations of methods for AffineSchemeOpenSubscheme                              #
########################################################################

########################################################################
# Intersections                                                        #
########################################################################
function intersect(
    Y::AbsAffineScheme,
    U::AffineSchemeOpenSubscheme;
    check::Bool=true
  )
  X = ambient_scheme(U)
  ambient_coordinate_ring(U) === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  X === Y && return AffineSchemeOpenSubscheme(Y, complement_equations(U), check=check)
  if check && !is_subscheme(Y, X)
    Y = intersect(Y, X)
  end
  return AffineSchemeOpenSubscheme(Y, [g for g in complement_equations(U) if !iszero(OO(Y)(g))], check=check)
end

intersect(U::AffineSchemeOpenSubscheme, Y::AbsAffineScheme) = intersect(Y, U)

function intersect(
    U::AffineSchemeOpenSubscheme,
    V::AffineSchemeOpenSubscheme
  )
  X = ambient_scheme(U)
  X == ambient_scheme(V) || error("ambient schemes do not coincide")
  return AffineSchemeOpenSubscheme(X, [a*b for a in complement_equations(U) for b in complement_equations(V)])
end

########################################################################
# Unions                                                               #
########################################################################
function Base.union(U::AffineSchemeOpenSubscheme, V::AffineSchemeOpenSubscheme)
  ambient_scheme(U) == ambient_scheme(V) || error("the two open sets are not contained in the same ambient scheme")
  return AffineSchemeOpenSubscheme(ambient_scheme(U), vcat(complement_equations(U), complement_equations(V)))
end

########################################################################
# Containment and equality                                             #
########################################################################
function is_subscheme(
    Y::AbsAffineScheme,
    U::AffineSchemeOpenSubscheme
  )
  ambient_coordinate_ring(Y) === ambient_coordinate_ring(U) || return false
  is_subscheme(Y, ambient_scheme(U)) || return false
  return one(OO(Y)) in ideal(OO(Y), complement_equations(U))
end


function is_subscheme(
    U::AffineSchemeOpenSubscheme,
    Y::AbsAffineScheme
  )
  return all(is_subscheme(V, Y) for V in affine_patches(U))
end

function is_subscheme(U::AffineSchemeOpenSubscheme, V::AffineSchemeOpenSubscheme)
  ambient_coordinate_ring(U) === ambient_coordinate_ring(V) || return false
  Z = complement(V)
  # perform an implicit radical membership test (Rabinowitsch) that is way more
  # efficient than computing radicals.
  for g in complement_equations(U)
    isempty(hypersurface_complement(Z, g)) || return false
  end
  return true
  #return is_subscheme(complement(intersect(V, ambient_scheme(U))), complement(U))
end

# TODO: Where did the type declaration go?
function ==(U::AffineSchemeSubset, V::AffineSchemeSubset)
  ambient_coordinate_ring(U) === ambient_coordinate_ring(V) || return false
  return is_subscheme(U, V) && is_subscheme(V, U)
end

########################################################################
# Closures of AffineSchemeOpenSubschemes                                                #
########################################################################
@doc raw"""
    closure(U::AffineSchemeOpenSubscheme)

Compute the Zariski closure of an open set ``U ⊂ X``
where ``X`` is the affine ambient scheme of ``U``.
"""
function closure(U::AffineSchemeOpenSubscheme{<:StdAffineScheme})
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = saturated_ideal(modulus(OO(X)))
  I = saturation(I, ideal(R, complement_equations(U)))
  return subscheme(X, I)
end

function closure(U::AffineSchemeOpenSubscheme{AffineSchemeType}) where {AffineSchemeType<:AffineScheme{<:Ring, <:MPolyRing}}
  return ambient_scheme(U)
end

function closure(U::AffineSchemeOpenSubscheme{AffineSchemeType}) where {AffineSchemeType<:AffineScheme{<:Ring, <:MPolyQuoRing}}
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = modulus(OO(X))
  I = saturation(I, ideal(R, complement_equations(U)))
  return subscheme(X, I)
end

@doc raw"""
    closure(U::AffineSchemeOpenSubscheme, Y::AbsAffineScheme)

Compute the closure of ``U ⊂ Y``.
"""
function closure(
    U::AffineSchemeOpenSubscheme,
    Y::AbsAffineScheme; check::Bool=true
  )
  @check is_subscheme(U, Y) "the first set is not contained in the second"
  X = closure(U)
  return intersect(X, Y)
end

########################################################################
# Preimages of open sets under AffineSchemeMors                                #
########################################################################

function preimage(f::AbsAffineSchemeMor, V::AffineSchemeOpenSubscheme; check::Bool=true)
  @check is_subscheme(codomain(f), ambient_scheme(V)) "set is not guaranteed to be open in the codomain"
  new_gens = pullback(f).(complement_equations(V))
  return AffineSchemeOpenSubscheme(domain(f), lifted_numerator.(new_gens), check=check)
end


########################################################################
# Printing                                                             #
########################################################################

function Base.show(io::IO, ::MIME"text/plain", U::AffineSchemeOpenSubscheme)
  io = pretty(io)
  println(io, "Open subset")
  println(io, Indent(), "of ", Lowercase(), ambient_space(U))
  print(io, Dedent(), "complement to V(")
  join(io, gens(complement_ideal(U)), ", ")
  print(io, ")")
end

# For the printing of regular functions, we need details on the affine patches.
# In general, one could avoid those details by just stating what is its ambient
# space and complement (see printing above)
function _show_semi_compact(io::IO, U::AffineSchemeOpenSubscheme)
  io = pretty(io)
  println(io, "Open subset")
  c = ambient_coordinates(U)
  str = "["*join(c, ", ")*"]"
  print(io, Indent(), "of affine scheme with coordinate")
  length(c) != 1 && print(io, "s")
  println(io, " "*str)
  print(io, Dedent(), "complement to V(")
  join(io, gens(complement_ideal(U)), ", ")
  print(io, ")")
  if n_patches(U) > 0
    println(io)
    l = ndigits(n_patches(U))
    print(io, "covered by ", ItemQuantity(n_patches(U), "affine patch"))
    print(io, Indent())
    co_str = [""]
    for V in affine_patches(U)
      cV = ambient_coordinates(V)
      str = "["*join(cV, ", ")*"]"
      push!(co_str, str)
    end
    k = max(length.(co_str)...)
    for i in 1:n_patches(U)
      li = ndigits(i)
      V = affine_patches(U)[i]
      println(io)
      kV = length(co_str[i+1])
      print(io, " "^(l-li)*"$(i): "*co_str[i+1]*" "^(k-kV+3), Lowercase(), V)
    end
    print(io, Dedent())
  end
end

function Base.show(io::IO, U::AffineSchemeOpenSubscheme)
  show_coord = get(io, :show_coordinates, true)
  if get(io, :show_semi_compact, false)
    _show_semi_compact(io, U)
  elseif isdefined(U, :name)
    print(io, name(U))
  elseif is_terse(io)
    print(io, "Open subset of affine scheme")
  elseif get_attribute(U, :is_empty, false)
    io = pretty(io)
    print(io, "Empty open subset of ", Lowercase(), ambient_space(U))
  else
    io = pretty(io)
    print(io, "Complement to V(")
    print(io, join(gens(complement_ideal(U)), ", "), ")")
    if show_coord
      c = ambient_coordinates(U)
      str = "["*join(c, ", ")*"]"
      print(io, " in affine scheme with coordinate")
      length(c) != 1 && print(io, "s")
      print(io, " "*str)
    end
  end
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, U::AffineSchemeOpenSubscheme;
    ambient_map::AbsAffineSchemeMor=base_change(phi, ambient_scheme(U))[2] # the base change on the ambient scheme
  )
  Y = domain(ambient_map)
  pbf = pullback(ambient_map)
  h = pbf.(complement_equations(U))
  UU = AffineSchemeOpenSubscheme(Y, h)
  return UU, restrict(ambient_map, UU, U, check=true) # TODO: Set to false after testing
end

