
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
function Base.show(io::IO, X::AbsCoveredScheme)
  if has_name(X)
    print(io, name_of(X))
    return
  end
  print(io, "covered scheme with $(npatches(default_covering(X))) affine patches in its default covering")
end
