
########################################################################
# Implementations of methods for SpecOpen                              #
########################################################################

########################################################################
# Intersections                                                        #
########################################################################
function intersect(
    Y::AbsSpec, 
    U::SpecOpen;
    check::Bool=true
  )
  X = ambient_scheme(U)
  ambient_coordinate_ring(U) === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  X === Y && return SpecOpen(Y, complement_equations(U), check=check)
  if check && !issubset(Y, X)
    Y = intersect(Y, X)
  end
  return SpecOpen(Y, [g for g in complement_equations(U) if !iszero(OO(Y)(g))], check=check)
end

intersect(U::SpecOpen, Y::AbsSpec) = intersect(Y, U)

function intersect(
    U::SpecOpen,
    V::SpecOpen
  )
  X = ambient_scheme(U)
  X == ambient_scheme(V) || error("ambient schemes do not coincide")
  return SpecOpen(X, [a*b for a in complement_equations(U) for b in complement_equations(V)])
end

########################################################################
# Unions                                                               #
########################################################################
function Base.union(U::SpecOpen, V::SpecOpen)
  ambient_scheme(U) == ambient_scheme(V) || error("the two open sets are not contained in the same ambient scheme")
  return SpecOpen(ambient_scheme(U), vcat(complement_equations(U), complement_equations(V)))
end

########################################################################
# Containment and equality                                             #
########################################################################
function issubset(
    Y::AbsSpec,
    U::SpecOpen
  )
  ambient_coordinate_ring(Y) === ambient_coordinate_ring(U) || return false
  issubset(Y, ambient_scheme(U)) || return false
  return one(OO(Y)) in ideal(OO(Y), complement_equations(U))
end


function issubset(
    U::SpecOpen,
    Y::AbsSpec
  ) 
  return all(issubset(V, Y) for V in affine_patches(U))
end

function issubset(U::SpecOpen, V::SpecOpen)
  ambient_coordinate_ring(U) === ambient_coordinate_ring(V) || return false
  Z = complement(V)
  # perform an implicit radical membership test (Rabinowitsch) that is way more 
  # efficient than computing radicals.
  for g in complement_equations(U)
    isempty(hypersurface_complement(Z, g)) || return false
  end
  return true
  #return issubset(complement(intersect(V, ambient_scheme(U))), complement(U))
end

# TODO: Where did the type declaration go?
function ==(U::SpecSubset, V::SpecSubset)
  ambient_coordinate_ring(U) === ambient_coordinate_ring(V) || return false
  return issubset(U, V) && issubset(V, U)
end

########################################################################
# Closures of SpecOpens                                                #
########################################################################
@doc raw"""
    closure(U::SpecOpen)

Compute the Zariski closure of an open set ``U ⊂ X`` 
where ``X`` is the affine ambient scheme of ``U``.
"""
function closure(U::SpecOpen{<:StdSpec})
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = saturated_ideal(modulus(OO(X)))
  I = saturation(I, ideal(R, complement_equations(U)))
  return subscheme(X, I)
end

function closure(U::SpecOpen{SpecType}) where {SpecType<:Spec{<:Ring, <:MPolyRing}}
  return ambient_scheme(U)
end

function closure(U::SpecOpen{SpecType}) where {SpecType<:Spec{<:Ring, <:MPolyQuoRing}}
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = modulus(OO(X))
  I = saturation(I, ideal(R, complement_equations(U)))
  return subscheme(X, I)
end

@doc raw"""
    closure(U::SpecOpen, Y::AbsSpec)

Compute the closure of ``U ⊂ Y``.
"""
function closure(
    U::SpecOpen,
    Y::AbsSpec; check::Bool=true
  )
  @check issubset(U, Y) "the first set is not contained in the second"
  X = closure(U)
  return intersect(X, Y)
end

########################################################################
# Preimages of open sets under SpecMors                                #
########################################################################

function preimage(f::AbsSpecMor, V::SpecOpen; check::Bool=true)
  @check issubset(codomain(f), ambient_scheme(V)) "set is not guaranteed to be open in the codomain"
  new_gens = pullback(f).(complement_equations(V))
  return SpecOpen(domain(f), lifted_numerator.(new_gens), check=check)
end


########################################################################
# Printing                                                             #
########################################################################

function Base.show(io::IO, ::MIME"text/plain", U::SpecOpen)
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
function _show_semi_compact(io::IO, U::SpecOpen)
  io = pretty(io)
  println(io, "Open subset")
  c = ambient_coordinates(U)
  str = "["*join(c, ", ")*"]"
  print(io, Indent(), "of affine scheme with coordinate")
  length(c) > 1 && print(io, "s")
  println(io, " "*str)
  print(io, Dedent(), "complement to V(")
  join(io, gens(complement_ideal(U)), ", ")
  print(io, ")")
  if npatches(U) > 0
    println(io)
    l = ndigits(npatches(U))
    print(io, "covered by $(npatches(U)) affine patch")
    npatches(U) > 1 && print(io, "es")
    print(io, Indent())
    co_str = [""]
    for V in affine_patches(U)
      cV = ambient_coordinates(V)
      str = "["*join(cV, ", ")*"]"
      push!(co_str, str)
    end
    k = max(length.(co_str)...)
    for i in 1:npatches(U)
      li = ndigits(i)
      V = affine_patches(U)[i]
      println(io)
      kV = length(co_str[i+1])
      print(io, " "^(l-li)*"$(i): "*co_str[i+1]*" "^(k-kV+3), Lowercase(), V)
    end
    print(io, Dedent())
  end
end

function Base.show(io::IO, U::SpecOpen, show_coord::Bool = true)
  io = pretty(io)
  if isdefined(U, :name) 
    print(io, name(U))
  elseif get(io, :supercompact, false)
    print(io, "Open subset of affine scheme")
  elseif get_attribute(U, :is_empty, false)
    print(io, "Empty open subset of ", Lowercase(), ambient_space(U))
  else
    print(io, "Complement to V(")
    print(io, join(gens(complement_ideal(U)), ", "), ")")
    if show_coord
      c = ambient_coordinates(U)
      str = "["*join(c, ", ")*"]"
      print(io, " in affine scheme with coordinate")
      if length(c) > 1
        print(io, "s")
      end
      print(io, " "*str)
    end
  end
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, U::SpecOpen;
    ambient_map::AbsSpecMor=base_change(phi, ambient_scheme(U))[2] # the base change on the ambient scheme
  )
  Y = domain(ambient_map)
  pbf = pullback(ambient_map)
  h = pbf.(complement_equations(U))
  UU = SpecOpen(Y, h)
  return UU, restrict(ambient_map, UU, U, check=true) # TODO: Set to false after testing
end

