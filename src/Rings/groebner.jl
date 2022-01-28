export  groebner_basis, groebner_basis_with_transformation_matrix, leading_ideal, syzygy_generators

# groebner stuff #######################################################
@doc Markdown.doc"""
    groebner_assure(I::MPolyIdeal; complete_reduction::Bool = false)
    groebner_assure(I::MPolyIdeal, ord::MonomialOrdering; complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an ideal `I` in a multivariate polynomial ring this function assures that a
Groebner basis w.r.t. the given monomial ordering is attached to `I` in `I.gb`.
It *currently* also ensures that the basis is defined on the Singular side in
`I.gb.S`, but this should not be relied upon: use `singular_assure(I.gb)` before
accessing `I.gb.S`.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ, ["x","y"], ordering=:degrevlex)
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*y-3*x,y^3-2*x^2*y])
ideal(x*y - 3*x, -2*x^2*y + y^3)

julia> Oscar.groebner_assure(I)
3-element Vector{fmpq_mpoly}:
 x*y - 3*x
 y^3 - 6*x^2
 x^3 - 9//2*x

julia> I.gb
Oscar.BiPolyArray{fmpq_mpoly}(fmpq_mpoly[x*y - 3*x, y^3 - 6*x^2, x^3 - 9//2*x], Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y - 3*x, y^3 - 6*x^2, x^3 - 9//2*x), Multivariate Polynomial Ring in x, y over Rational Field, Singular Polynomial Ring (QQ),(x,y),(dp(2),C), true, #undef, false)
```
"""
function groebner_assure(I::MPolyIdeal; complete_reduction::Bool = false)
  if !isdefined(I, :gb)
    singular_assure(I)
    i = Singular.std(I.gens.S; complete_reduction = complete_reduction)
    I.gb = BiPolyArray(base_ring(I), i)
    I.gb.isGB = true
    I.gb.O = [I.gb.Ox(x) for x = gens(I.gb.S)]
  end
end

@doc Markdown.doc"""
    groebner_basis(B::BiPolyArray; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `BiPolyArray` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of the
ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`. The Groebner basis is then
returned in `B.S`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> A = Oscar.BiPolyArray([x*y-3*x,y^3-2*x^2*y])
Oscar.BiPolyArray{fmpq_mpoly}(fmpq_mpoly[x*y - 3*x, -2*x^2*y + y^3], #undef, Multivariate Polynomial Ring in x, y over Rational Field, #undef, false, #undef, true)

julia> B = groebner_basis(A)
Oscar.BiPolyArray{fmpq_mpoly}(fmpq_mpoly[#undef, #undef, #undef], Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y - 3*x, y^3 - 6*x^2, 2*x^3 - 9*x), Multivariate Polynomial Ring in x, y over Rational Field, Singular Polynomial Ring (QQ),(x,y),(dp(2),C), true, #undef, true)
```
"""
function groebner_basis(B::BiPolyArray; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
  # if ord != :degrevlex
    R = singular_ring(B.Ox, ordering)
    i = Singular.Ideal(R, [R(x) for x = B])
#    @show "std on", i, B
    i = Singular.std(i, complete_reduction = complete_reduction)
    return BiPolyArray(B.Ox, i)
  # end
  if !isdefined(B, :S)
    B.S = Singular.Ideal(B.Sx, [B.Sx(x) for x = B.O])
  end
#  @show "dtd", B.S
  return BiPolyArray(B.Ox, Singular.std(B.S, complete_reduction = complete_reduction))
end

#= @doc Markdown.doc"""
 =     groebner_basis(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
 =
 = Compute a Groebner basis w.r.t. the given monomial ordering of the polynomial ring.
 =
 = # Examples
 = ```jldoctest
 = julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
 = (Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])
 =
 = julia> I = ideal([x*y-3*x,y^3-2*x^2*y])
 = ideal(x*y - 3*x, -2*x^2*y + y^3)
 =
 = julia> G = groebner_basis(I)
 = 3-element Vector{fmpq_mpoly}:
 =  x*y - 3*x
 =  y^3 - 6*x^2
 =  x^3 - 9//2*x
 = ```
 = """ =#
#= function groebner_basis(I::MPolyIdeal)
 =   groebner_assure(I)
 =   return collect(I.gb)
 = end =#

@doc Markdown.doc"""
    groebner_basis(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
    groebner_basis(I::MPolyIdeal, ord::MonomialOrdering; complete_reduction::Bool=false)
 
Given an ideal `I` and optional parameters monomial ordering `ordering` and `complete_reduction`,
compute a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of `I`
    w.r.t. the given monomial ordering `ordering` (as default `:degrevlex`).

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*y-3*x,y^3-2*x^2*y])
ideal(x*y - 3*x, -2*x^2*y + y^3)

julia> H = groebner_basis(I; ordering=:lex)
3-element Vector{fmpq_mpoly}:
 y^4 - 3*y^3
 x*y - 3*x
 -y^3 + 6*x^2
```
"""
function groebner_basis(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
  R = singular_ring(base_ring(I), ordering)
  !Oscar.Singular.has_global_ordering(R) && error("The ordering has to be a global ordering.")
  i = Singular.std(Singular.Ideal(R, [R(x) for x = gens(I)]), complete_reduction = complete_reduction)
  return collect(BiPolyArray(base_ring(I), i))
end

@doc Markdown.doc"""
    groebner_basis_with_transform(B::BiPolyArray; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
    groebner_basis_with_transform(B::BiPolyArray, ord::MonomialOrdering; complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `BiPolyArray` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of the
ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering` and the transformation matrix from the ideal to the Groebner basis. Return value is a BiPolyArray together with a map.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> A = Oscar.BiPolyArray([x*y-3*x,y^3-2*x^2*y])
Oscar.BiPolyArray{fmpq_mpoly}(fmpq_mpoly[x*y - 3*x, -2*x^2*y + y^3], #undef, Multivariate Polynomial Ring in x, y over Rational Field, #undef, false, #undef, true)

julia> B,m = Oscar.groebner_basis_with_transform(A)
(Oscar.BiPolyArray{fmpq_mpoly}(fmpq_mpoly[#undef, #undef, #undef], Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y - 3*x, y^3 - 6*x^2, 6*x^3 - 27*x), Multivariate Polynomial Ring in x, y over Rational Field, Singular Polynomial Ring (QQ),(x,y),(dp(2),C), false, #undef, true), [1 2*x -2*x^2+y^2+3*y+9; 0 1 -x])
```
"""
function groebner_basis_with_transform(B::BiPolyArray; ord::Symbol = :degrevlex, complete_reduction::Bool = false)
  if ord != :degrevlex
    R = singular_ring(B.Ox, ord)
    i = Singular.Ideal(R, [R(x) for x = B])
#    @show "std on", i, B
    i, m = Singular.lift_std(i, complete_reduction = complete_reduction)
    return BiPolyArray(B.Ox, i), map_entries(x->B.Ox(x), m)
  end
  if !isdefined(B, :S)
    singular_assure(B)
  end
#  @show "dtd", B.S

  i, m = Singular.lift_std(B.S, complete_reduction = complete_reduction)
  return BiPolyArray(B.Ox, i), map_entries(x->B.Ox(x), m)
end

@doc Markdown.doc"""
    groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool=false)
    groebner_basis_with_transformation_matrix(I::MPolyIdeal, ord::MonomialOrdering; complete_reduction::Bool=false)

Return a pair `G, m` where `G` is a Groebner basis of the ideal `I` with respect to the
monomial ordering `ordering`, and `m` is a transformation matrix from `gens(I)` to `G`. If
`complete_reduction` is set to `true` then `G` will be the reduced Groebner basis.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ,["x","y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*y^2-1,x^3+y^2+x*y])
ideal(x*y^2 - 1, x^3 + x*y + y^2)

julia> G,m = groebner_basis_with_transformation_matrix(I)
(fmpq_mpoly[x*y^2 - 1, x^3 + x*y + y^2, x^2 + y^4 + y], fmpq_mpoly[1 0; 0 1; -x^2 - y y^2])

julia> m * gens(I) == G
true
```
"""
function groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool=false)
  G, m = Oscar.groebner_basis_with_transform(I; ordering=ordering, complete_reduction=complete_reduction)
  return G, Array(m)
end

# syzygies #######################################################
@doc Markdown.doc"""
    syzygy_generators(a::Vector{<:MPolyElem})

Return generators for the syzygies on the given polynomials.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> S = syzygy_generators([x^3+y+2,x*y^2-13*x^2,y-14])
3-element Vector{FreeModElem{fmpq_mpoly}}:
 (-y + 14)*e[2] + (-13*x^2 + x*y^2)*e[3]
 (-169*y + 2366)*e[1] + (-13*x*y + 182*x - 196*y + 2744)*e[2] + (13*x^2*y^2 - 2548*x^2 + 196*x*y^2 + 169*y + 338)*e[3]
 (-13*x^2 + 196*x)*e[1] + (-x^3 - 16)*e[2] + (x^4*y + 14*x^4 + 13*x^2 + 16*x*y + 28*x)*e[3]
```
"""
function syzygy_generators(a::Vector{<:MPolyElem})
  I = ideal(a)
  singular_assure(I)
  s = Singular.syz(I.gens.S)
  F = free_module(parent(a[1]), length(a))
  @assert rank(s) == length(a)
  return [F(s[i]) for i=1:Singular.ngens(s)]
end

# leading ideal #######################################################
@doc Markdown.doc"""
    leading_ideal(g::Vector{T}, args...) where { T <: MPolyElem }

Return the ideal generated by the leading monomials of the given polynomials.
If not otherwise given as a further argument this is done w.r.t. the
degree reverse lexicographical monomial ordering.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5])
ideal(x*y^2, x^3)
```
"""
function leading_ideal(g::Vector{T}, args...) where { T <: MPolyElem }
  return ideal([ leading_monomial(f, args...) for f in g ])
end

@doc Markdown.doc"""
    leading_ideal(g::Vector{Any}, args...)

Given an array of fitting type this function returns the ideal
generated by the leading monomials of the given array elements.
If not otherwise given as a further argument this is done w.r.t.
the degree reverse lexicographical monomial ordering.
"""
function leading_ideal(g::Vector{Any}, args...)
  return leading_ideal(typeof(g[1])[ f for f in g ], args...)
end

@doc Markdown.doc"""
    leading_ideal(Rx::MPolyRing, g::Vector{Any}, args...)

Given a multivariate polynomial ring `Rx` and an array of elements
this function generates an array of multivariate polynomials in `Rx`
and returns the leading ideal for the ideal generated by the given
array elements w.r.t. the given monomial ordering of `Rx`.

"""
function leading_ideal(Rx::MPolyRing, g::Vector{Any}, args...)
  h = elem_type(Rx)[ Rx(f) for f in g ]
  return leading_ideal(h, args...)
end

@doc Markdown.doc"""
    leading_ideal(I::MPolyIdeal)
    leading_ideal(I::MPolyIdeal, ord::MonomialOrdering)

Given a multivariate polynomial ideal `Ì` this function returns the
leading ideal for `I`. This is done w.r.t. the given monomial ordering
in the polynomial ring of `I`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R,[x*y^2-3*x, x^3-14*y^5])
ideal(x*y^2 - 3*x, x^3 - 14*y^5)

julia> L = leading_ideal(I)
ideal(x*y^2, x^4, y^5)
```
"""
function leading_ideal(I::MPolyIdeal)
  singular_assure(I)
  groebner_assure(I)
  singular_assure(I.gb)
  return MPolyIdeal(base_ring(I), Singular.Ideal(I.gb.Sx, [Singular.leading_monomial(g) for g in gens(I.gb.S)]))
end
 
function leading_ideal(I::MPolyIdeal, ord::MonomialOrdering)
  singular_assure(I, ord)
  groebner_assure(I, ord)
  singular_assure(I.gb, ord)
  return MPolyIdeal(base_ring(I), Singular.Ideal(I.gb.Sx, [Singular.leading_monomial(g) for g in gens(I.gb.S)]))
end
 
@doc Markdown.doc"""
    leading_ideal(I::MPolyIdeal, ordering::Symbol)

Given a multivariate polynomial ideal `Ì` and a monomial ordering `ordering`
this function returns the leading ideal for `I` w.r.t. `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R,[x*y^2-3*x, x^3-14*y^5])
ideal(x*y^2 - 3*x, x^3 - 14*y^5)

julia> L = leading_ideal(I, :lex)
ideal(y^7, x*y^2, x^3)
```
"""
function leading_ideal(I::MPolyIdeal, ordering::Symbol)
  return leading_ideal(groebner_basis(I; ordering=ordering), ordering)
end

@doc Markdown.doc"""
    normal_form_internal(I::Singular.sideal, J::MPolyIdeal)

**Note**: Internal function, subject to change, do not use.

Compute the normal form of the generators `gens(I)` of the ideal `I` w.r.t. a
Groebner basis of `J`.

CAVEAT: This computation needs a Groebner basis of `J`. If this Groebner basis
is not available, one is computed automatically. This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> groebner_basis(J)
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> SR = singular_ring(base_ring(J))
Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C)

julia> I = Singular.Ideal(SR,[SR(-1+c+b+a^3),SR(-1+b+c*a+2*a^3),SR(5+c*b+c^2*a)])
Singular ideal over Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C) with generators (a^3 + b + c - 1, 2*a^3 + a*c + b - 1, a*c^2 + b*c + 5)

julia> Oscar.normal_form_internal(I,J)
3-element Vector{fmpq_mpoly}:
 a^3
 2*a^3 + 2*a - 2*c
 4*a - 2*c^2 - c + 5
```
"""
function normal_form_internal(I::Singular.sideal, J::MPolyIdeal)
    if !isdefined(J, :gb)
        groebner_assure(J)
    end
    K = ideal(base_ring(J), reduce(I, J.gb.S))
    return [J.gens.Ox(x) for x = gens(K.gens.S)]
end

@doc Markdown.doc"""
    normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }

Compute the normal form of the polynomial `f` w.r.t. a
Groebner basis of `J`.

CAVEAT: This computation needs a Groebner basis of `J`. If this Groebner basis
is not available, one is computed automatically. This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> groebner_basis(J)
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(-1+c+b+a^3, J)
a^3
```
"""
function normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }
    singular_assure(J)
    I = Singular.Ideal(J.gens.Sx, J.gens.Sx(f))
    N = normal_form_internal(I, J)
    return N[1]
end

@doc Markdown.doc"""
    normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }

Compute the normal form of the elements of the array `A` w.r.t. a
Groebner basis of `J`.

CAVEAT: This computation needs a Groebner basis of `J`. If this Groebner basis
is not available, one is computed automatically. This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> A = [-1+c+b+a^3,-1+b+c*a+2*a^3,5+c*b+c^2*a]
3-element Vector{fmpq_mpoly}:
 a^3 + b + c - 1
 2*a^3 + a*c + b - 1
 a*c^2 + b*c + 5

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> groebner_basis(J)
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(A, J)
3-element Vector{fmpq_mpoly}:
 a^3
 2*a^3 + 2*a - 2*c
 4*a - 2*c^2 - c + 5
```
"""
function normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }
    singular_assure(J)
    I = Singular.Ideal(J.gens.Sx, [J.gens.Sx(x) for x in A])
    normal_form_internal(I, J)
end

function groebner_assure(I::MPolyIdeal, ord::MonomialOrdering; complete_reduction::Bool = false)
  R = base_ring(I)
  Rx = singular_ring(R, ord.o)

  if !isdefined(I, :gb) || ordering(I.gb.Sx) != ordering(Rx)
    I.gens.Sx = Rx
    I.gens.S = Singular.Ideal(Rx, Rx.(gens(I)))
    I.gb = BiPolyArray(I.gens.Ox, Singular.std(I.gens.S, complete_reduction = complete_reduction))
  end
end

function groebner_basis(B::BiPolyArray, ord::MonomialOrdering; complete_reduction::Bool = false)
   singular_assure(B, ord)
   R = B.Sx
   !Oscar.Singular.has_global_ordering(R) && error("The ordering has to be a global ordering.")
   I = Singular.Ideal(R, gens(B.S)...)
   i = Singular.std(I, complete_reduction = complete_reduction)
   return BiPolyArray(B.Ox, i)
end

function groebner_basis(I::MPolyIdeal, ord::MonomialOrdering; complete_reduction::Bool=false)
  return collect(groebner_basis(I.gens, ord))
end

function groebner_basis_with_transform(B::BiPolyArray, ord::MonomialOrdering; complete_reduction::Bool = false)
   if !isdefined(B, :ord)
      singular_assure(B, ord)
   elseif ord != B.ord
     R = singular_ring(B.Ox, ord)
     i = Singular.Ideal(R, [R(x) for x = B])
     i, m = Singular.lift_std(i, complete_reduction = complete_reduction)
     return BiPolyArray(B.Ox, i), map_entries(x->B.Ox(x), m)
   end

   if !isdefined(B, :S)
     B.S = Singular.Ideal(B.Sx, [B.Sx(x) for x = B.O])
   end
 
   i, m = Singular.lift_std(B.S, complete_reduction = complete_reduction)
   return BiPolyArray(B.Ox, i), map_entries(x->B.Ox(x), m)
 end
 
 function groebner_basis_with_transformation_matrix(I::MPolyIdeal, ord::MonomialOrdering; complete_reduction::Bool=false)
   G, m = Oscar.groebner_basis_with_transform(I, ord; complete_reduction=complete_reduction)
   return G, Array(m)
 end
