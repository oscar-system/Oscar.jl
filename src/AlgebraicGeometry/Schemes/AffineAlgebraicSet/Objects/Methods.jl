########################################################
# (1) Display
########################################################

# As a set, corresponding the fat ideal or the radical does not change anything.
# If one already knows the reduced scheme structure on this algebraic set, we
# print it (more convenient to see the ideal)
function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Affine algebraic set")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  if isdefined(X, :Xred)
    I = saturated_ideal(defining_ideal(X))
  else
    I = fat_ideal(X)
  end
  print(io, Dedent(), "defined by ", Lowercase(), I)
end

# As a set, corresponding the fat ideal or the radical does not change anything.
# If one already knows the reduced scheme structure on this algebraic set, we
# print it (more convenient to see the ideal)
#
# For compact printing, we value the notation V(bla) since it tells everything
# we need to know, in a given contextual printing
function Base.show(io::IO, X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  if is_terse(io)
    print(io, "Affine algebraic set")
  elseif get_attribute(X, :is_empty, false)
    io = pretty(io)
    print(io, "Empty affine algebraic set over ")
    K = base_ring(X)
    print(terse(io), Lowercase(), K)
  else
    io = pretty(io)
    if isdefined(X, :Xred)
      I = saturated_ideal(defining_ideal(X))
    else
      I = fat_ideal(X)
    end
    print(io, LowercaseOff(), "V(")
    join(io, gens(I), ", ")
    print(io,")")
  end
end

# special case for Zariski opens
function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet)
  io = pretty(io)
  println(io, "Reduced subscheme")
  print(io, Indent(),"of ", Lowercase(), fat_scheme(X))
  print(io, Dedent())
end

function Base.show(io::IO, X::AffineAlgebraicSet)
  if is_terse(io)
    print(io, "Affine algebraic set")
  elseif get_attribute(X, :is_empty, false)
    io = pretty(io)
    print(io, "Empty affine algebraic set over ")
    K = base_ring(X)
    print(terse(io), Lowercase(), K)
  else
    io = pretty(io)
    print(io, "Reduced subscheme of ", Lowercase(), fat_scheme(X))
  end
end

