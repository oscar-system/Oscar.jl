########################################################################
# Check for emptyness                                                  #
########################################################################

@attr Bool function Base.isempty(U::AffineSchemeOpenSubscheme)
  return all(isempty, affine_patches(U))
end

########################################################################
# Methods for VarietyFunctionField                                     #
########################################################################
### essential getters
representative_patch(KK::VarietyFunctionField) = KK.U
variety(KK::VarietyFunctionField) = KK.X
# provide one alias for the user's convenience
scheme(KK::VarietyFunctionField) = variety(KK)
coefficient_ring(KK::VarietyFunctionField) = KK.kk
representative_field(KK::VarietyFunctionField) = KK.KK

### user facing constructors


@doc raw"""
    function_field(X::AbsCoveredScheme)

Return the function field of the irreducible variety `X`.

Internally, a rational function is represented by an element in the field of
fractions of the `ambient_coordinate_ring` of the `representative_patch`.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P, I);

julia> Ycov = covered_scheme(Y)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme(-(y//x)^2*(z//x) + 1)
    2: scheme((x//y)^3 - (z//y))
    3: scheme((x//z)^3 - (y//z)^2)
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> K = function_field(Ycov)
Field of rational functions
  on scheme over QQ covered with 3 patches
    1: [(y//x), (z//x)]   scheme(-(y//x)^2*(z//x) + 1)
    2: [(x//y), (z//y)]   scheme((x//y)^3 - (z//y))
    3: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2)
represented by
  patch 1: fraction field of multivariate polynomial ring

julia> one(K)
Rational function
  on scheme over QQ covered with 3 patches
    1: [(y//x), (z//x)]   scheme(-(y//x)^2*(z//x) + 1)
    2: [(x//y), (z//y)]   scheme((x//y)^3 - (z//y))
    3: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2)
represented by
  patch 1: 1
```
"""
@attr VarietyFunctionField function function_field(X::AbsCoveredScheme; check::Bool=true)
  return VarietyFunctionField(X, check=check)
end

########################################################################
# Methods for VarietyFunctionFieldElem                                 #
########################################################################
### essential getters 
numerator(f::VarietyFunctionFieldElem) = numerator(f.f)
denominator(f::VarietyFunctionFieldElem) = denominator(f.f)
parent(f::VarietyFunctionFieldElem) = f.KK
representative(f::VarietyFunctionFieldElem) = f.f

### constructors
one(KK::VarietyFunctionField) = VarietyFunctionFieldElem(KK, one(ambient_coordinate_ring(representative_patch(KK))),
                                                         one(ambient_coordinate_ring(representative_patch(KK)))
                                                        )
zero(KK::VarietyFunctionField) = VarietyFunctionFieldElem(KK, zero(ambient_coordinate_ring(representative_patch(KK))),
                                                          one(ambient_coordinate_ring(representative_patch(KK)))
                                                         )

### arithmetic 
function +(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(representative(a) + representative(b), check=false)
end

function -(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(representative(a) - representative(b), check=false)
end

function -(a::T) where {T<:VarietyFunctionFieldElem}
  return (parent(a))((-1)*representative(a), check=false)
end

function *(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(representative(a) * representative(b), check=false)
end

function Base.:(//)(a::Integer, b::T) where {T<:VarietyFunctionFieldElem}
  return (parent(b))(a//representative(b))
end

function Base.:(//)(a::ZZRingElem, b::T) where {T<:VarietyFunctionFieldElem}
  return (parent(b))(a//representative(b))
end

function Base.:(//)(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  return (parent(a))(representative(a) // representative(b), check=false)
end

function ==(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  KK = parent(a)
  U = representative_patch(KK)
  return iszero(OO(U)(numerator(a)*denominator(b) - numerator(b)*denominator(a)))
end

# We need to manually split this into three methods, because 
# otherwise it seems that Julia can not dispatch this function.
function ^(a::VarietyFunctionFieldElem, i::Int64)
  return parent(a)(representative(a)^i, check=false)
end
function ^(a::VarietyFunctionFieldElem, i::Integer)
  return parent(a)(representative(a)^i, check=false)
end
function ^(a::VarietyFunctionFieldElem, i::ZZRingElem)
  return parent(a)(representative(a)^i, check=false)
end

# Multiplication with elements of the coefficient ring.
# Eventually this needs to be restricted further in case 
# of ambiguities.
function *(a::RingElem, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === coefficient_ring(parent(b)) || error("the arguments do not have the same parent ring")
  return (parent(b))(a * representative(b), check=false)
end
function *(b::T, a::RingElem) where {T<:VarietyFunctionFieldElem}
  return a*b
end

function *(a::Integer, b::T) where {T<:VarietyFunctionFieldElem}
  return coefficient_ring(parent(b))(a)*b
end

function *(b::T, a::Integer) where {T<:VarietyFunctionFieldElem}
  return a*b
end

# try to avoid a groebner basis computation
iszero(a::VarietyFunctionFieldElem) = iszero(representative(a)) || iszero(OO(representative_patch(parent(a)))(numerator(a)))
isone(a::VarietyFunctionFieldElem) = isone(representative(a)) || iszero(OO(representative_patch(parent(a)))(numerator(a) - denominator(a)))
is_unit(a::VarietyFunctionFieldElem) = !iszero(representative(a))

########################################################################
# Conversion of rational functions on arbitrary patches                #
########################################################################

function (KK::VarietyFunctionField)(a::RingElem; check::Bool=true)
  return KK(a, one(parent(a)), check=check)
end



function (KK::VarietyFunctionField)(a::MPolyQuoRingElem, b::MPolyQuoRingElem; check::Bool=true)
  return KK(lift(a), lift(b), check=check)
end

function (KK::VarietyFunctionField)(a::MPolyQuoLocRingElem, 
                                    b::MPolyQuoLocRingElem; 
                                    check::Bool=true)
  return KK(lifted_numerator(a)*lifted_denominator(b), lifted_denominator(a)*lifted_numerator(b), check=check)
end

function (KK::VarietyFunctionField)(a::MPolyLocRingElem, 
                                    b::MPolyLocRingElem; 
                                    check::Bool=true)
  return KK(numerator(a)*denominator(b), denominator(a)*numerator(b), check=check)
end

function (KK::VarietyFunctionField)(a::MPolyRingElem, b::MPolyRingElem; check::Bool=true)
  R = parent(a)
  R === parent(b) || error("rings are not compatible")
  R === ambient_coordinate_ring(representative_patch(KK)) && return VarietyFunctionFieldElem(KK, a, b)
  
  # otherwise check whether we can find the ring of h among the affine patches
  if R in [ambient_coordinate_ring(V) for V in patches(default_covering(variety(KK)))] 
    # allocate a variable for the patch in which a and be are living
    V = representative_patch(KK)
    X = variety(KK)
    C = default_covering(X)
    for i in 1:n_patches(C)
      if ambient_coordinate_ring(C[i]) === R
        V = C[i]
        break
      end
    end

    # convert it 
    U = representative_patch(KK)
    h_generic = move_representative(a, b, V, U, C)
    return VarietyFunctionFieldElem(KK, numerator(h_generic),
                                    denominator(h_generic)
                                   )
  else
    # go through the registered coverings and look for the ring
    X = variety(KK)
    for C in coverings(X)
      for U in patches(C)
        if R === ambient_coordinate_ring(U)
          iso = _flatten_open_subscheme(U, default_covering(X))
          pb = pullback(inverse(iso))
          pb_num = fraction(pb(OO(U)(a)))
          pb_den = fraction(pb(OO(U)(b)))
          h = pb_num//pb_den # a representative on an affine chart of X
          V = representative_patch(KK)
          W = ambient_scheme(codomain(iso))
          if ambient_coordinate_ring(V) === base_ring(parent(h))
            return VarietyFunctionFieldElem(KK, numerator(h), denominator(h))
          else
            hh = move_representative(numerator(h), denominator(h), W, V, default_covering(X))
            return VarietyFunctionFieldElem(KK, numerator(hh), denominator(hh))
          end
        end
      end
    end
  end

  error("no open subset found with admissible ring")
end

@doc raw"""
    function move_representative(
        a::MPolyRingElem, b::MPolyRingElem,
        V::AbsAffineScheme, U::AbsAffineScheme,
        C::Covering
      )

Given a fraction ``a/b ∈ Quot(P)`` with ``P = 𝕜[x]`` the `ambient_coordinate_ring`
of an affine patch ``V`` in a covering `C`, move that fraction to 
one in ``Quot(P')`` where ``P'`` is the ambient coordinate ring of another patch ``U``.

**Note:** This is only guaranteed to work for irreducible schemes! 
"""
function move_representative(
    a::MPolyRingElem, b::MPolyRingElem,
    V::AbsAffineScheme, U::AbsAffineScheme,
    C::Covering
  )
  G = C[U, V]
  f, _ = gluing_morphisms(G)
  A, B = gluing_domains(G)
  pba = pullback(f)(OO(B)(a))
  pbb = pullback(f)(OO(B)(b))
  iszero(pbb) && error("pullback of denominator is zero")
  # in the next line, A is either a AffineSchemeOpenSubscheme or a PrincipalOpenSubset
  h_generic = generic_fraction(pba, A)//generic_fraction(pbb, A)
  return h_generic

  # It turned out that the following is too expensive in general
  if domain(f) isa PrincipalOpenSubset
    fac = factor(lifted_numerator(complement_equation(domain(f))))
    p = OO(U)(numerator(h_generic))
    q = OO(U)(denominator(h_generic))
    for (a, e) in fac
      aa = OO(U)(a)
      k_num, _ = _minimal_power_such_that(aa, x->divides(p, x)[1]) # This division takes ages for big polynomials.
      k_den, _ = _minimal_power_such_that(aa, x->divides(q, x)[1])
      k = minimum([k_num, k_den])
      aa = aa^k
      _, p = divides(p, aa)
      _, q = divides(q, aa)
    end
    h_generic = fraction(p)//fraction(q)
  end
  return h_generic
end

function (KK::VarietyFunctionField)(h::AbstractAlgebra.Generic.FracFieldElem; check::Bool=true)
  return KK(numerator(h), denominator(h), check=check)
end

### given the fraction field of the `ambient_coordinate_ring` in one
# affine chart, return a representative of `f` in that field
function (K::AbstractAlgebra.Generic.FracField)(f::VarietyFunctionFieldElem)
  R = base_ring(K)
  V = representative_patch(parent(f))
  C = default_covering(variety(parent(f)))
  for U in patches(C)
    if ambient_coordinate_ring(U) == R
      V = U
      break
    end
  end
  f_mov = move_representative(numerator(f), denominator(f), 
                              representative_patch(parent(f)),
                              V, C
                             )
  return K(numerator(f_mov), denominator(f_mov))
end

function getindex(f::VarietyFunctionFieldElem, V::AbsAffineScheme)
  C = default_covering(variety(parent(f))) 
  if any(x->x===V, patches(C))
    return move_representative(numerator(f), denominator(f), 
                               representative_patch(parent(f)),
                               V, C
                              )
  else
    KK = parent(f)
    X = variety(KK)
    iso = _flatten_open_subscheme(V, C)
    VV = codomain(iso)
    W = ambient_scheme(VV)
    g = f[W] # Will go to the first case above
    num = pullback(iso)(OO(VV)(numerator(g)))
    den = pullback(iso)(OO(VV)(denominator(g)))
    return fraction(num)//fraction(den)
  end
end

# some dummy methods for compatibility
fraction(a::MPolyQuoRingElem) = lift(a)//one(lift(a))
fraction(a::MPolyRingElem) = a//one(a)

########################################################################
# Implementation of the rest of the interfaces                         #
########################################################################

function elem_type(::Type{T}) where {FracFieldType, 
                                     T<:VarietyFunctionField{<:Field, 
                                                             <:FracFieldType
                                                            }
                                    }
  return VarietyFunctionFieldElem{elem_type(FracFieldType), T}
end

function parent_type(::Type{T}) where {ParentType, T<:VarietyFunctionFieldElem{<:Any, ParentType}}
  return ParentType
end

base_ring(KK::VarietyFunctionField) = KK.kk
base_ring_type(::Type{T}) where {R, T<:VarietyFunctionField{R}} = R

is_domain_type(::Type{T}) where {T<:VarietyFunctionFieldElem} = true
is_exact_type(::Type{T}) where {T<:VarietyFunctionFieldElem} = true

function Base.hash(f::VarietyFunctionFieldElem, h::UInt)
  r = 57103
  return xor(r, hash(representative(f), h))
end

function Base.deepcopy_internal(f::VarietyFunctionFieldElem, dict::IdDict)
  return parent(f)(deepcopy_internal(numerator(representative(f)), dict), 
                   deepcopy_internal(denominator(representative(f)), dict), 
                   check=false)
end

(KK::VarietyFunctionField)() = zero(KK)
function (KK::VarietyFunctionField)(a::Integer)
  R = base_ring(representative_field(KK))
  KK(R(a), one(R), check=false)
end
(KK::VarietyFunctionField)(f::VarietyFunctionFieldElem) = (parent(f) == KK ? f : error("element does not belong to the given field"))
(KK::VarietyFunctionField)(a::MPolyRingElem) = KK(a, one(parent(a)), check=false)
canonical_unit(f::VarietyFunctionFieldElem) = f # part of the ring interface that becomes trivial for fields

function Base.show(io::IO, KK::VarietyFunctionField)
  io = pretty(io)
  if is_terse(io)
    print(io, "Field of rational functions")
  else
    print(io, "Field of rational functions on ", Lowercase(), variety(KK))
  end
end

# The function field is global, but we know how it is represented on a given
# chart of the covering of X: once we have described and labeled the charts, we
# mention on which chart we have that representation
function Base.show(io::IO, ::MIME"text/plain", KK::VarietyFunctionField)
  io = pretty(io)
  X = variety(KK)
  cov = default_covering(X)
  println(io, "Field of rational functions")
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => cov), X)
  j = findfirst(U -> representative_patch(KK) === U, collect(cov))
  println(io, Dedent())
  println(io, "represented by")
  print(io, Indent(), "patch $j: ", Lowercase(), representative_field(KK))
  print(io, Dedent())
end

function Base.show(io::IO, f::VarietyFunctionFieldElem)
  io = pretty(io)
  X = variety(parent(f))
  if get(io, :show_semi_compact, false)
    cov = Oscar._covering_for_printing(io, X)
    k = get(io, :offset, 0)
    _show_semi_compact(io, f, cov, k)
  elseif is_terse(io)
    print(io, "Rational function")
  else
    print(io, "Rational function on ", Lowercase(), variety(parent(f)))
  end
end

# Needed in nested printings where we already know the parent, but need only a
# description of the rational function on the appropriate charts
function _show_semi_compact(io::IO, f::VarietyFunctionFieldElem, cov::Covering, k::Int)
  io = pretty(io)
  j = findfirst(U -> representative_patch(parent(f)) === U, collect(cov))
  print(io, "Rational function represented by ", representative(f), " "^k, " on patch $j")
end

# Same details as for the printing of the parent
function Base.show(io::IO, ::MIME"text/plain", f::VarietyFunctionFieldElem)
  io = pretty(io)
  KK = parent(f)
  X = variety(KK)
  cov = default_covering(X)
  j = findfirst(U -> representative_patch(KK) === U, collect(cov))
  println(io, "Rational function")
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => cov), X)
  println(io, Dedent())
  println(io, "represented by")
  print(io, Indent(), "patch $j: ", representative(f))
  print(io, Dedent())
end

function divexact(f::VarietyFunctionFieldElem, 
    g::VarietyFunctionFieldElem;
    check::Bool=true
  )
  return f//g
end
inv(f::VarietyFunctionFieldElem) = parent(f)(denominator(representative(f)),
                                      numerator(representative(f)),
                                      check=false
                                     )

AbstractAlgebra.promote_rule(::Type{T}, ::Type{S}) where {T<:VarietyFunctionFieldElem, S<:Integer} = T
AbstractAlgebra.promote_rule(::Type{T}, ::Type{S}) where {T<:VarietyFunctionFieldElem, S<:AbstractAlgebra.Generic.FracFieldElem} = T

AbstractAlgebra.promote_rule(::Type{T}, ::Type{T}) where {T<:VarietyFunctionFieldElem} = T

function AbstractAlgebra.promote_rule(::Type{FFET}, ::Type{U}) where {T, FFET<:VarietyFunctionFieldElem{T}, U<:RingElement}
  promote_rule(T, U) == T ? FFET : Union{}
end 

function AbstractAlgebra.promote_rule(::Type{FFET}, ::Type{U}) where {T, FFET<:VarietyFunctionFieldElem{AbstractAlgebra.Generic.FracFieldElem{T}}, U<:RingElement}
  promote_rule(T, U) == T ? FFET : Union{}
end 

characteristic(KK::VarietyFunctionField) = characteristic(base_ring(variety(KK)))

@doc raw"""
    is_regular(f::VarietyFunctionFieldElem, U::Scheme)

Return whether ``f ∈ K(X)`` restricts to a regular function 
on an open subset ``U ⊂ X``.

**Note:** ``U`` must either be an `affine_chart` of ``X`` or 
its `ambient_scheme` must be an `affine_chart`.
"""
function is_regular(f::VarietyFunctionFieldElem, U::Scheme)
  error("method not implemented")
end

function is_regular(f::VarietyFunctionFieldElem, U::AbsAffineScheme)
  p = numerator(f[U])
  q = denominator(f[U])
  return _is_regular_fraction(OO(U), p, q)
end

function is_regular(f::VarietyFunctionFieldElem, U::PrincipalOpenSubset)
  Y = ambient_scheme(U)
  ff = f[Y]
  p = numerator(ff)
  q = denominator(ff)
  return _is_regular_fraction(OO(U), p, q)
end

function is_regular(f::VarietyFunctionFieldElem, W::AffineSchemeOpenSubscheme)
  return all(U->is_regular(f, U), affine_patches(W))
end

function pullback(f::AbsCoveredSchemeMorphism, a::VarietyFunctionFieldElem)
  fcov = covering_morphism(f)
  KK = parent(a)
  V = representative_patch(KK)
  if any(U->codomain(fcov[U])===V, patches(domain(fcov)))
    j = findfirst(U->codomain(fcov[U])===V, patches(domain(fcov)))
    U = patches(domain(fcov))[j]
    phi = pullback(fcov[U])
    return function_field(domain(f))(phi(OO(V)(numerator(a))), 
                                     phi(OO(V)(denominator(a))))
  else 
    U = first(patches(domain(fcov)))
    W = codomain(fcov[U])
    #(any(x->x===V, affine_charts(codomain(f))) && any(x->x===W, affine_charts(codomain(f)))) || error("method not implemented for too complicated transitions")
    b = a[W]
    phi = pullback(fcov[U])
    return function_field(domain(f))(phi(OO(W)(numerator(b))), 
                                     phi(OO(W)(denominator(b))))
  end
end

scheme(f::VarietyFunctionFieldElem) = scheme(parent(f))
variety(f::VarietyFunctionFieldElem) = variety(parent(f))

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(K::VarietyFunctionField)
  F = representative_field(K)
  P = base_ring(F)::MPolyRing
  num = rand(P, 0:5, 0:5, 0:5)
  den = zero(P)
  while is_zero(den)
    den = rand(P, 1:5, 1:5, 1:5)
  end
  return K(num, den)
end
