export set_name!


####################################################################################
# (1) Equality
####################################################################################

# TODO: Add further cross-type comparison methods as needed.

function ==(X::T, Y::T) where {T<:AbsSpec{<:Ring, <:MPolyRing}}
  return OO(X) === OO(Y)
end


function ==(X::AbsSpec, Y::AbsSpec)
  X === Y && return true
  return issubset(X, Y) && issubset(Y, X)
end


function ==(X::AbsSpec, Y::EmptyScheme)
  return issubset(X, Y)
end


==(X::EmptyScheme, Y::AbsSpec) = (Y == X)



########################################################
# (2) Display
########################################################

function Base.show(io::IO, X::AbsSpec)
  if has_attribute(X, :name)
    print(io, name(X))
    return
  end
  print(io, "Spec of $(OO(X))")
end



########################################################
# (3) Check for zero divisors in rings
########################################################

@Markdown.doc """
    is_non_zero_divisor(f::RingElem, X::AbsSpec)

Checks if a ring element is a non-zero divisor
in the coordinate ring of an affine scheme.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1, x2, x3) = gens(OO(X))
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> is_non_zero_divisor(x1, X)
true

julia> is_non_zero_divisor(zero(OO(X)), X)
false
```
"""
function is_non_zero_divisor(f::RingElem, X::AbsSpec)
  error("method not implemented for affine schemes of type $(typeof(X))")
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::MPolyQuoElem, X::AbsSpec{<:Ring, <:MPolyQuo})
  R = ambient_coordinate_ring(X)
  I = modulus(OO(X))
  J = ideal(R, lift(f))
  return I == quotient(I, J)
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = ideal(OO(X), [zero(OO(X))])
  zero_ideal = Oscar.pre_image_ideal(I)
  J = Oscar.pre_image_ideal(ideal(OO(X), [f]))
  Q = quotient(zero_ideal, J)
  return zero_ideal == Q
end

########################################################################
# High level constructors of subschemes                                #
########################################################################

function sub(X::AbsSpec, I::Ideal)
  inc = ClosedEmbedding(X, I) # Will complain if input is not compatible
  return domain(inc), inc
end

########################################################################
# Further intersections                                                #
########################################################################

function intersect(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
  if !(ambient_scheme(U) === ambient_scheme(V))
    return intersect(underlying_scheme(U), underlying_scheme(V))
  end
  return PrincipalOpenSubset(ambient_scheme(U), complement_equation(U)*complement_equation(V))
end


# TODO: The method below is necessary for a temporary hotfix; see #1882.
# It seems likely that the implementation can be tuned so that, for 
# instance, the massive use of primary decompositions can be avoided. 
# This should eventually be addressed. 
@Markdown.doc """
    components(X::AbsSpec)

This returns a decomposition of ``X`` into its connected components 
``X = U₁ ∪ U₂ ∪ … ∪ Uₙ`` with ``Uᵢ`` a `PrincipalOpenSubset` of ``X``.
"""
function components(X::AbsSpec)
  I = saturated_ideal(modulus(OO(X)))
  l = primary_decomposition(I)
  comp = [subscheme(X, Q) for (Q, _) in l]
  found_intersection = true
  while found_intersection
    found_intersection = false
    for i in 1:length(comp)-1
      for j in i+1:length(comp)
        if !isempty(intersect(comp[i], comp[j]))
          found_intersection = true
          C2 = popat!(comp, j)
          C1 = comp[i]
          C_new = subscheme(X, intersect(saturated_ideal(modulus(OO(C1))),
                                         saturated_ideal(modulus(OO(C2)))
                                        ))
          comp[i] = C_new
        end
      end
    end
  end

  # We need to reproduce the components of X as `PrincipalOpenSubset`s. 
  # To this end, we first determine separating equations between pairs 
  # of components.

  n = length(comp)
  separating_equations = Array{RingElem}(undef, n, n)
  for i in 1:length(comp)
    C1 = comp[i]
    for j in i:length(comp)
      C2 = comp[j]
      if i == j
        separating_equations[i, j] = one(OO(X))
        continue
      end
      
      v = OO(X).(gens(saturated_ideal(modulus(OO(C1)))))
      w = OO(X).(gens(saturated_ideal(modulus(OO(C2)))))

      J = ideal(OO(X), vcat(v, w))
      # The idea is to write 1 = f + g with f ∈ I(C1) and g ∈ I(C2).
      # Then 1 - f = g vanishes identically on C2, but is a unit in OO(C1) 
      # and vice versa. 
      c = coordinates(one(OO(X)), J)

      # Conversion to assure compatibility. 
      # This can be removed, once the ideal interface is streamlined. 
      if c isa MatrixElem
        c = [c[1, i] for i in 1:ncols(c)]::Vector
      end

      a = c[1:length(v)]
      b = c[length(v)+1:length(v)+length(w)]

      f = sum([v[i]*a[i] for i in 1:length(v)], init=zero(OO(X)))
      g = sum([w[i]*b[i] for i in 1:length(w)], init=zero(OO(X)))
      separating_equations[i, j] = 1 - f
      separating_equations[j, i] = 1 - g
    end
  end
  result = Vector{PrincipalOpenSubset}()
  for i in 1:n
    f = prod(separating_equations[i, :], init=one(OO(X)))
    push!(result, PrincipalOpenSubset(X, f))
  end
  return result
end




