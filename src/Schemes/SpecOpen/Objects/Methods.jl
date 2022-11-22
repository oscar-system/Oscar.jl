export intersect

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
  X === Y && return SpecOpen(Y, gens(U), check=check)
  if check && !issubset(Y, X)
    Y = intersect(Y, X)
  end
  return SpecOpen(Y, gens(U), check=check)
end

intersect(U::SpecOpen, Y::AbsSpec) = intersect(Y, U)

function intersect(
    U::SpecOpen,
    V::SpecOpen
  )
  X = ambient_scheme(U)
  X == ambient_scheme(V) || error("ambient schemes do not coincide")
  return SpecOpen(X, [a*b for a in gens(U) for b in gens(V)])
end

########################################################################
# Unions                                                               #
########################################################################
function Base.union(U::SpecOpen, V::SpecOpen)
  ambient_scheme(U) == ambient_scheme(V) || error("the two open sets are not contained in the same ambient scheme")
  return SpecOpen(ambient_scheme(U), vcat(gens(U), gens(V)))
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
  return one(OO(Y)) in ideal(OO(Y), gens(U))
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
  for g in gens(U)
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
@Markdown.doc """
    closure(U::SpecOpen)

Compute the Zariski closure of an open set ``U ⊂ X`` 
where ``X`` is the affine ambient scheme of ``U``.
"""
function closure(U::SpecOpen{<:StdSpec})
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = saturated_ideal(modulus(OO(X)))
  I = saturation(I, ideal(R, gens(U)))
  return subscheme(X, I)
end

function closure(U::SpecOpen{SpecType}) where {SpecType<:Spec{<:Ring, <:MPolyRing}}
  return ambient_scheme(U)
end

function closure(U::SpecOpen{SpecType}) where {SpecType<:Spec{<:Ring, <:MPolyQuo}}
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = modulus(OO(X))
  I = saturation(I, ideal(R, gens(U)))
  return subscheme(X, I)
end

@Markdown.doc """
    closure(U::SpecOpen, Y::AbsSpec)

Compute the closure of ``U ⊂ Y``.
"""
function closure(
    U::SpecOpen,
    Y::AbsSpec 
  )
  issubset(U, Y) || error("the first set is not contained in the second")
  X = closure(U)
  return intersect(X, Y)
end

########################################################################
# Preimages of open sets under SpecMors                                #
########################################################################

function preimage(f::AbsSpecMor, V::SpecOpen; check::Bool=true)
  if check
    issubset(codomain(f), ambient_scheme(V)) || error("set is not guaranteed to be open in the codomain")
  end
  new_gens = pullback(f).(gens(V))
  return SpecOpen(domain(f), lifted_numerator.(new_gens), check=check)
end


########################################################################
# Printing                                                             #
########################################################################
function Base.show(io::IO, U::SpecOpen)
  if isdefined(U, :name) 
    print(io, name(U))
    return
  end
  print(io, "complement of zero locus of $(gens(U)) in $(ambient_scheme(U))")
end

