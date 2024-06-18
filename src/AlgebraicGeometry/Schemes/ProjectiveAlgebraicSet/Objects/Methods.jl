########################################################
# (1) Display
########################################################

function Base.show(io::IO, ::MIME"text/plain",
                   X::AbsProjectiveAlgebraicSet{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Projective algebraic set")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  if isdefined(X, :Xred)
    I = vanishing_ideal(X)
  else
    I = fat_ideal(X)
  end
  print(io, Dedent(), "defined by ", Lowercase(), I)
end

# If we know a radical ideal describing our algebraic set, we preferably print
# that one (it is easier to read...)
function Base.show(io::IO, X::AbsProjectiveAlgebraicSet{<:Field, <:MPolyQuoRing})
  if is_terse(io)
    print(io, "Projective algebraic set")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty projective algebraic set")
  else
    io = pretty(io)
    print(io, LowercaseOff(), "V(")
    if isdefined(X, :Xred)
      I = vanishing_ideal(X)
    else
      I = fat_ideal(X)
    end
    join(io, gens(I), ",")
    print(io, ")")
  end
end

