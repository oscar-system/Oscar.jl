
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
  if length(C) > 0
    println(io, Dedent(), "with default covering")
    print(io, Indent(), "described by patch")
    print(io, Indent())
      for i in 1:length(C)
      println(io)
      print(io, "$i: ", Lowercase(), C[i])
    end
    println(io)
    print(io, Dedent(), "in the coordinate(s)")
    print(io, Indent())
    for i in 1:length(C)
      println(io)
      print(io, "$i: [")
      co = coordinates(C[i])
      if length(co) == 0
        print(io, "]")
        continue
      end
      for j in 1:length(co)-1
        print(io, "$(co[j]), ")
      end
      print(io, "$(co[end])]")
    end
    print(io, Dedent())
    print(io, Dedent())
  else
    print(io, "with empty default covering")
  end
end

function _show_semi_compact(io::IO, X::AbsCoveredScheme, cov::Covering = get_attribute(X, :simplified_covering, default_covering(X)), n::String = "")
  io = pretty(io)
  show(io, X, cov)
  print(io, Indent())
  co_str = String[]
  for i in 1:length(cov)
    U = cov[i]
    co = coordinates(U)
    str = reduce(*, ["$x, " for x in co], init = "[")
    str = str[1:end-2]*"]"
    push!(co_str, str)
  end
  k = max(length.(co_str)...)
  for i in 1:length(cov)
    U = cov[i]
    kc = length(co_str[i])
    println(io)
    print(io, "$(i)"*n*": "*co_str[i]*" "^(k-kc+3), Lowercase(), U)
  end
  print(io, Dedent())
end

function Base.show(io::IO, X::AbsCoveredScheme, cov::Covering = get_attribute(X, :simplified_covering, default_covering(X)))
  io = pretty(io)
  n = npatches(cov)
  if has_name(X)
    print(io, name(X))
  elseif get(io, :supercompact, false)
    print(io, "Scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty covered scheme")
  else
    print(io, "Scheme over ")
    if base_ring(X) == QQ
      print(io, "QQ")
    else
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
    end
    print(io, " covered with $n patch")
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
