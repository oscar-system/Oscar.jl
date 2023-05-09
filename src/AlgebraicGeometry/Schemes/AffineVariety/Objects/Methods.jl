########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyQuoRing})
  println(io, "Affine variety")
  print(io, " in ")
  println(io, ambient_space(X))
  print(io, "defined by ")
  print(io, vanishing_ideal(X))
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyQuoRing})
# one line printing
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  else
    print(io, "Affine variety defined by $(vanishing_ideal(X))")
  end
end


# The case that X is not Zariski closed in its ambient space.
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety)
  println(io, "Affine variety")
  print(io, "  in ")
  println(io, ambient_space(X))
  print(io, "with coordinate ring $(OO(X))")
end

function Base.show(io::IO, X::AffineVariety)
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  else
    print(io, "Affine variety contained in the vanishing locus of $(vanishing_ideal(X))")
  end
end


# For affine space
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyRing})
  println(io, "Affine space of dimension $(dim(X))")
  print(io, "  with coordinates ")
  for x in coordinates(X)
    print(io, x, " ")
  end
  println(io, "")
  print(io, "  over ")
  print(io, base_ring(X))
end


function Base.show(io::IO, X::AffineVariety{<:Field, <:MPolyRing})
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
      print(IOContext(io, :supercompact => true), base_ring(X))
    else
      # nested printing allowed, preferably supercompact
      print(io, "Affine $(dim(X))-space over ")
      print(IOContext(io, :supercompact => true), base_ring(X))
    end
  end
end
