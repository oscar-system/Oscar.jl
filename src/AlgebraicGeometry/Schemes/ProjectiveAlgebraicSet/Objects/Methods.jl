########################################################
# (1) Display
########################################################

function Base.show(io::IO, ::MIME"text/plain",
                   X::AbsProjectiveAlgebraicSet{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Algebraic set")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  print(io, Dedent(), "defined by ", Lowercase(), defining_ideal(X))
end

# If we know a radical ideal describing our algebraic set, we preferably print
# that one (it is easier to read...)
function Base.show(io::IO, X::AbsProjectiveAlgebraicSet{<:Field, <:MPolyQuoRing})
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty projective algebraic set")
  else
    print(io, LowercaseOff(), "V(")
    if isdefined(X, :Xred)
      I = vanishing_ideal(X)
    else
      I = fat_ideal(X)
    end
    join(io, gens(I), ",")
    print(io,")")
  end
end

