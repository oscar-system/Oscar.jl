
########################################################################
# Required methods for AbsCoveredScheme                                #
########################################################################
### Forwarding of the non-documented getters
Base.in(U::AbsSpec, X::AbsCoveredScheme) = (U in underlying_scheme(X))::Bool
getindex(X::AbsCoveredScheme, i::Int) = coverings(underlying_scheme(X))[i]::Covering
getindex(X::AbsCoveredScheme, C::Covering, D::Covering) = getindex(underlying_scheme(X), C, D)::CoveringMorphism
#setindex(X::AbsCoveredScheme, f::CoveringMorphism, C::Covering, D::Covering) = setindex(underlying_scheme(X), f, C, D)::CoveringMorphism
glueings(X::AbsCoveredScheme) = glueings(underlying_scheme(X))::IdDict{<:Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing}

########################################################################
# Finding coverings, patches, and glueings                             #
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

Base.in(U::AbsSpec, X::CoveredScheme) = any(C->(U in C), coverings(X))

########################################################################
# Printing                                                             #
########################################################################

function Base.show(io::IO, ::MIME"text/plain", X::AbsCoveredScheme)
  io = pretty(io)
  println(io, "Scheme")
  println(io, Indent(), "over ", Lowercase(), base_ring(X))
  C = default_covering(X)
  n = npatches(C)
  print(io, Dedent(), "with $n affine patch")
  n > 1 && print(io, "es")
  print(io, " in its default covering")
  print(io, Indent())
  for U in C
    println(io)
    print(io, Lowercase(), U)
  end
  print(io, Dedent())
end

function _show_semi_compact(io::IO, X::AbsCoveredScheme, cov::Covering = get_attribute(X, :simplified_covering, default_covering(X)), l::Int = 0)
  io = pretty(io)
  show(io, X, cov)
  n = npatches(cov)
  print(io, Indent())
  for U in cov
    println(io)
    print(io, " "^l, Lowercase(), U)
  end
  print(io, Dedent())
end

function Base.show(io::IO, X::AbsCoveredScheme, cov::Covering = get_attribute(X, :simplified_covering, default_covering(X)))
  io = pretty(io)
  n = npatches(cov)
  if has_name(X)
    print(io, name(X))
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty covered scheme")
  elseif get(io, :supercompact, false)
    print(io, "Scheme")
  else
    print(io, "Scheme over ")
    if base_ring(X) == QQ
      print(io, "QQ")
    else
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
    end
    print(io, " covered with $n affine patch")
    n > 1 && print(io, "es")
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
