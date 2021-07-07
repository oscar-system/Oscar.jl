##############################################################################
#
# Mapping between rings preserving variable names. The behaviour is essentially
# that of printing the polynomial in the domain, then evaluating the string
# in the codomain. Of course this is supposed to be done efficiently.
#
##############################################################################

export embed_names

##############################################################################
#
# 1. implementations for R1[x...] -> R2[x...]
#
##############################################################################

mutable struct IMapMPolyToMPoly{S, T, SE}
  D::S
  R::T
  var_map::Vector{Int}
  have_missing::Bool    # purely for optimization purposes
end

domain(m::IMapMPolyToMPoly) = m.D
codomain(m::IMapMPolyToMPoly) = m.R

function (m::IMapMPolyToMPoly{S, T, SE})(d::SE) where
                              {S <: MPolyRing, T <: MPolyRing, SE <: MPolyElem}
  parent(d) == m.D || error("parent mismatch")
  r = MPolyBuildCtx(m.R)
  nD = ngens(m.D)
  nR = ngens(m.R)
  re = zeros(Int, nR)
  cR = coefficient_ring(m.R)
  have_missing = m.have_missing
  var_map = m.var_map
  for (dc, de) in zip(coefficients(d), exponent_vectors(d))
    for j in 1:nR
      re[j] = 0
    end
    if have_missing
      for i in 1:nD
        de[i] > 0 || continue
        if var_map[i] < 0
          # variable i appears in this term of d and is mapped to zero
          @goto skip
        end
        re[var_map[i]] += de[i]
      end
    else
      for i in 1:nD
        re[var_map[i]] += de[i]
      end
    end
    push_term!(r, cR(dc), re)
@label skip
  end
  return finish(r)
end

function (m::IMapMPolyToMPoly{S, T, SE})(d::MPolyIdeal{SE}) where
                              {S <: MPolyRing, T <: MPolyRing, SE <: MPolyElem}
  return ideal(m.R, [m(i) for i in gens(d)])
end


@doc Markdown.doc"""
    embed_names(D::S, C::T) where {S <: MPolyRing, T <: MPolyRing}

Return a function that maps elements of $D$ to elements of $C$ by the names
of the variables. If a variable in $D$ does not appear in $C$, it maps to zero.
The target $C$ must not have duplicate variables names among the names
of the variables in $D$.
"""
function embed_names(D::S, R::T) where {S <: MPolyRing, T <: MPolyRing}
  sD = symbols(D)
  sR = symbols(R)
  var_map = [-1 for i in 1:length(sD)]
  have_missing = false
  for i in 1:length(sD)
    inds = findall(x -> x == sD[i], sR)
    if length(inds) == 1
      var_map[i] = inds[1]
    elseif length(inds) > 1
      error("duplicate variable $(sD[i]) found among $sR")
    else
      have_missing = true
    end
  end
  return IMapMPolyToMPoly{S, T, elem_type(S)}(D, R, var_map, have_missing)
end


##############################################################################
#
# 2. implementations for Q[t..., x...] -> Q(t...)[x...]
#
##############################################################################


