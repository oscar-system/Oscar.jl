########################################################
# (1) Display
########################################################

# detailed printing: if the variety if reduced, we preferably print the
# a radical ideal.. Otherwise, we do not compute anything and we print an ideal
# whose radical is the vanishing of X
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Affine variety")
  print(io, Indent(), "in ")
  println(io, Lowercase(), ambient_space(X))
  if get_attribute(X, :is_reduced, false)
    I = saturated_ideal(defining_ideal(X))
  else
    I = fat_ideal(X)
  end
  print(io, Dedent(), "defined by ", Lowercase(), I)
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyQuoRing})
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  elseif get_attribute(X, :is_empty, false)
    io = pretty(io)
    print(io, "Empty affine variety over ")
    K = base_ring(X)
    print(IOContext(io, :supercompact => true), Lowercase(), K)
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
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  elseif get_attribute(X, :is_empty, false)
    io = pretty(io)
    print(io, "Empty affine variety over ")
    K = base_ring(X)
    print(IOContext(io, :supercompact => true), Lowercase(), K)
  else
    io = pretty(io)
    print(io, "Affine open subset of ", Lowercase(), closure(X, ambient_space(X)))
  end
end


# For affine space: we print the coordinates, can be useful
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyRing})
  io = pretty(io)
  println(io, "Affine space of dimension $(dim(X))")
  print(io, Indent(), "over ")
  println(io, Lowercase(), base_ring(X))
  print(io, Dedent(), "with coordinate")
  length(coordinates(X)) != 1 && print(io, "s")
  print(io, " [")
  print(io, join(coordinates(X), ", "), "]")
end

# In a more compact printing, we allow unicode printing, whenever unicode is
# allowed, but we keep as an option of not printing the coordinates. For
# instance, in some nested printings for morphisms, we might have already
# mentioned the coordinates, and so `show_coord = false` ensures that we avoid
# redundancy
function Base.show(io::IO, X::AffineVariety{<:Field, <:MPolyRing})
  io = pretty(io)
  show_coord = get(io, :show_coordinates, true)
  if get(io, :supercompact, false)
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, LowercaseOff(), "𝔸$(ltx["\\^$(dim(X))"])")
    else
      print(io, LowercaseOff(), "AA^$(dim(X))")
    end
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine space over ")
    K = base_ring(X)
    print(IOContext(io, :supercompact => true), Lowercase(), K)
  else
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, LowercaseOff(), "𝔸")
      n = dim(X)
      for d in reverse(digits(n))
        print(io, ltx["\\^$d"])
      end
      print(io, " over ")
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
      if show_coord
        c = coordinates(X)
        print(io, " with coordinate")
        length(c) != 1 && print(io, "s")
        print(io, " [")
        print(io, join(c, ", "), "]")
      end
    else
      print(io, "Affine $(dim(X))-space over ")
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
      if show_coord
        c = coordinates(X)
        print(io, " with coordinate")
        length(c) != 1 && print(io, "s")
        print(io, " [")
        print(io, join(c, ", "), "]")
      end
    end
  end
end
