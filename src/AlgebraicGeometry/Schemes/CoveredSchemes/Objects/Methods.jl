
########################################################################
# Required methods for AbsCoveredScheme                                #
########################################################################
### Forwarding of the non-documented getters
Base.in(U::AbsAffineScheme, X::AbsCoveredScheme) = (U in underlying_scheme(X))::Bool
getindex(X::AbsCoveredScheme, i::Int) = coverings(underlying_scheme(X))[i]::Covering
getindex(X::AbsCoveredScheme, C::Covering, D::Covering) = getindex(underlying_scheme(X), C, D)::CoveringMorphism
#setindex(X::AbsCoveredScheme, f::CoveringMorphism, C::Covering, D::Covering) = setindex(underlying_scheme(X), f, C, D)::CoveringMorphism
gluings(X::AbsCoveredScheme) = gluings(underlying_scheme(X))::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:AbsGluing}

########################################################################
# Finding coverings, patches, and gluings                             #
########################################################################
getindex(X::CoveredScheme, C::CoveringType, D::CoveringType) where {CoveringType<:Covering} = X.refinements[(C, D)]
setindex(X::CoveredScheme, f::CoveringMorphismType, C::CoveringType, D::CoveringType) where {CoveringMorphismType<:CoveringMorphism, CoveringType<:Covering} = X.refinements[(C, D)]
getindex(X::CoveredScheme, i::Int) = coverings(X)[i]


function getindex(X::AbsCoveredScheme, C::Covering)
  for i in 1:length(coverings(X))
    C == coverings(X)[i] && return i
  end
  error("covering not listed")
end

getindex(X::AbsCoveredScheme, i::Int, j::Int) = X[X[i], X[j]]

Base.in(U::AbsAffineScheme, X::CoveredScheme) = any(C->(U in C), coverings(X))

########################################################################
# Printing                                                             #
########################################################################

# We use several trickes here to ensure a nice alignment of the information.
# First, we agreed on printing the detailed on the covering of the scheme - one
# may always be interested in the local information about the (covered) scheme.
#
# We label each patch so that in the case of many (>= 5) patches, it is easier
# to see to what patch we refer (and the ordering match with the ordering of
# `default_covering(X)` here).
#
# We first then create a list of strings for all the coordinates, for each
# chart, and look for the length of the longest string: this allows us to
# estimate an offset to add everything we print the coordinates. Hence, in front
# of each coordinate system, one prints the description of the charts and with a
# good handling of the offsets, the description should all be aligned on the
# lefts, 3 spaces after the longest set of coordinates.
#
# If we have more than 10 patches, we need to take care of a left offset, so
# that the label of the patches and aligned on the right: we use here the
# variables `l` and `li` which check the number of digits
function Base.show(io::IO, ::MIME"text/plain", X::AbsCoveredScheme)
  io = pretty(io)
  println(io, "Scheme")
  println(io, Indent(), "over ", Lowercase(), base_ring(X))
  C = default_covering(X)
  if length(C) > 0
    l = ndigits(length(C))
    println(io, Dedent(), "with default covering")
    print(io, Indent(), "described by patches")
    print(io, Indent())
      for i in 1:length(C)
      li = ndigits(i)
      println(io)
      print(io, " "^(l-li)*"$i: ", Lowercase(), C[i])
    end
    println(io)
    print(io, Dedent(), "in the coordinate(s)")
    print(io, Indent())
    for i in 1:length(C)
      li = ndigits(i)
      println(io)
      print(io, " "^(l-li)*"$i: ")
      co = coordinates(C[i])
      str = "["*join(co, ", ")*"]"
      print(io, str)
    end
    print(io, Dedent())
    print(io, Dedent())
  else
    print(io, "with empty default covering")
  end
end

_covering_for_printing(io::IO, X::AbsCoveredScheme) = get(io, :covering, get_attribute(X, :simplified_covering, default_covering(X)))

# The explanation of the printing are the same as the function above. We need
# this `_show_semi_compact` for detailed nested printing. Here, we assume that
# in our nest, we already made precised the base ring of our covered scheme.
# 
# We then only print detailed on the covering. Note that in some cases, we night
# not want to only print the default covering, but a simplified one or a very
# specific one (when one prints ideal sheaves, for instance). In that case, we
# need the flexibility on the covering, using the variable `cov`.
#
# When one has to deal with morphisms of covered schemes, we want to distinguish
# charts of the domain and chart of the codomain. So, in addition to the
# numbering of the labels, we add the possibility to have more using the string
# `n`. Note that, in default use in CoveredSchemeMor printing, `n` is either "a"
# or "b" depending on whether `X` is the domain or the codomain.
function _show_semi_compact(io::IO, X::AbsCoveredScheme, cov::Covering, n::String)
  io = pretty(io)
  show(IOContext(io, :covering => cov, :show_semi_compact => false), X)
  print(io, Indent())
  co_str = String[""]
  for i in 1:length(cov)
    U = cov[i]
    co = coordinates(U)
    str = "["*join(co, ", ")*"]"
    push!(co_str, str)
  end
  k = max(length.(co_str)...)
  for i in 1:length(cov)
    l = ndigits(length(cov))
    li = ndigits(i)
    U = cov[i]
    kc = length(co_str[i+1])
    println(io)
    print(io, " "^(l-li)*"$(i)"*n*": "*co_str[i+1]*" "^(k-kc+3), Lowercase(), U)
  end
  print(io, Dedent())
end

function Base.show(io::IO, X::AbsCoveredScheme)
  cov = _covering_for_printing(io, X)
  if get(io, :show_semi_compact, false)
    l = get(io, :label, "")
    _show_semi_compact(io, X, cov, l)
  else
    io = pretty(io)
    n = n_patches(cov)
    if has_name(X)
      print(io, name(X))
    elseif get(io, :supercompact, false)
      print(io, "Covered scheme")
    else
      if get_attribute(X, :is_empty, false) || n == 0
        print(io, "Empty covered scheme over ")
      else
        print(io, "Scheme over ")
      end
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
    end
    n > 0 && print(io, " covered with ", ItemQuantity(n, "patch"))
  end
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, X::AbsCoveredScheme)
  C = default_covering(X)
  CC, f_CC = base_change(phi, C)
  XX = CoveredScheme(CC)
  return XX, CoveredSchemeMorphism(XX, X, f_CC)
end
