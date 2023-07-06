########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Affine variety")
  print(io, Indent(), "in ")
  println(io, Lowercase(), ambient_space(X))
  print(io, "defined by ")
  print(io, Lowercase(), fat_ideal(X))
  print(io, Dedent())
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyQuoRing})
# one line printing
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  else
    print(io, X.X)
  end
end


# The case that X is not Zariski closed in its ambient space.
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety)
  io = pretty(io)
  println(io, "Affine variety")
  print(io, Indent(), "in ")
  println(io, ambient_space(X))
  print(io, "with coordinate ring ", Lowercase(), OO(X))
  print(io, Dedent())
end

function Base.show(io::IO, X::AffineVariety)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  else
    print(io, "Affine open in ", Lowercase(), closure(X))
  end
end


# For affine space
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyRing})
  io = pretty(io)
  println(io, "Affine space of dimension $(dim(X))")
  print(io, Indent(), "with coordinates ")
  for x in coordinates(X)
    print(io, x, " ")
  end
  println(io, "")
  print(io, "over ")
  print(io, Lowercase(), base_ring(X))
  print(io, Dedent())
end


function Base.show(io::IO, X::AffineVariety{<:Field, <:MPolyRing})
  io = pretty(io)
  if get(io, :supercompact, false) # no nested printing
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, "𝔸$(ltx["\\^$(dim(X))"])")
    else
      print(io, "AA^$(dim(X))")
    end
  else
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, "𝔸")
      n = dim(X)
      for d in reverse(digits(n))
        print(io, ltx["\\^$d"])
      end
      print(io, " over ")
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
    else
      # nested printing allowed, preferably supercompact
      print(io, "Affine $(dim(X))-space over ")
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
    end
  end
end
