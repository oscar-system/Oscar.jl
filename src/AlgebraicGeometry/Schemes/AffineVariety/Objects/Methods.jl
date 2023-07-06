########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Affine variety")
  print(io, Indent(), "in ")
  println(io, Lowercase(), ambient_space(X))
  print(io, Dedent(), "defined by ")
  print(io, Lowercase(), fat_ideal(X))
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine variety")
  else
    print(io, underlying_scheme(X))
  end
end


# The case that X is not Zariski closed in its ambient space.
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety)
  io = pretty(io)
  println(io, "Affine variety")
  print(io, Indent(), "in ")
  println(io, ambient_space(X))
  print(io, Dedent(), "with coordinate ring ", Lowercase(), OO(X))
end

function Base.show(io::IO, X::AffineVariety)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine variety")
  else
    print(io, "Affine open subset of ", Lowercase(), closure(X))
  end
end


# For affine space
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyRing})
  io = pretty(io)
  println(io, "Affine space of dimension $(dim(X))")
  print(io, Indent(), "over ")
  println(io, Lowercase(), base_ring(X))
  print(io, Dedent(), "with coordinate")
  length(coordinates(X)) > 1 && print(io, "s")
  print(io, " ")
  print(io, join(coordinates(X), ", "))
end


function Base.show(io::IO, X::AffineVariety{<:Field, <:MPolyRing})
  io = pretty(io)
  if get(io, :supercompact, false)
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, "𝔸$(ltx["\\^$(dim(X))"])")
    else
      print(io, "AA^$(dim(X))")
    end
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine space")
  else
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, "𝔸")
      n = dim(X)
      for d in reverse(digits(n))
        print(io, ltx["\\^$d"])
      end
      print(io, " over ")
      if base_ring(X) == QQ
        print(io, "QQ")
      else
        print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
      end
      c = coordinates(X)
      print(io, " with coordinate")
      length(c) > 1 && print(io, "s")
      print(io, " ")
      print(io, join(c, ", "))
    else
      print(io, "Affine $(dim(X))-space over ")
      if base_ring(X) == QQ
        print(io, "QQ")
      else
        print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
      end
      c = coordinates(X)
      print(io, " with coordinate")
      length(c) > 1 && print(io, "s")
      print(io, " ")
      print(io, join(c, ", "))
    end
  end
end
