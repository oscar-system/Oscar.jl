# TODO: in this file are used many methods TEMPORARILY defined in files matrix_manipulation.jl and stuff_field_gen.jl
# once methods in those files will be deleted / replaced / modified, this file need to be modified too



export
    generalized_jordan_block,
    generalized_jordan_form,
    is_conjugate,
    is_semisimple,
    is_unipotent,
    pol_elementary_divisors,
    multiplicative_jordan_decomposition




########################################################################
#
# Semisimple / Unipotent
#
########################################################################

"""
    multiplicative_jordan_decomposition(M::MatrixGroupElem)

Return `S` and `U` in the group `G = parent(M)` such that `S` is semisimple,
`U` is unipotent and  `M = SU = US`.
!!! warning "WARNING:" 
    this is *NOT*, in general, the same output returned when `M` has type `MatElem`.
"""
function multiplicative_jordan_decomposition(x::MatrixGroupElem)
   a = order(x)
   p = characteristic(base_ring(x))
   alpha = valuation(a,p)
   m = div(a, p^alpha)
   k = crt(fmpz(0),fmpz(p^alpha),fmpz(1),fmpz(m))
#   a,b = multiplicative_jordan_decomposition(x.elm)
   return x^k, x^(a+1-k)
end

"""
    is_semisimple(x::MatrixGroupElem{T}) where T <: FinFieldElem

Return whether `x` is semisimple, i.e. has order coprime with the characteristic of its base ring.
"""
is_semisimple(x::MatrixGroupElem{T}) where T <: FinFieldElem = is_coprime(order(Int, x), Int(characteristic(x.parent.ring)))

"""
    is_unipotent(x::MatrixGroupElem{T}) where T <: FinFieldElem

Return whether `x` is unipotent, i.e. its order is a power of the characteristic of its base ring.
"""
is_unipotent(x::MatrixGroupElem{T}) where T <: FinFieldElem = isone(x) || is_power(order(Int, x))[2]==Int(characteristic(x.parent.ring))






########################################################################
#
# MyIsConjugate
#
########################################################################



# return an element in the centralizer of x in GL(n,F) with determinant d
# Method: compute generators for the centralizer of x in GL(n,F);
# then, multiply them in order to get an element of determinant d.
# x has type MatrixGroupElem
function _elem_given_det(x,d)
   C,e = centralizer(GL(x.parent.deg, x.parent.ring),x)
   U,fa = unit_group(x.parent.ring)
   GA,ea = sub(U, [preimage(fa,det(g)) for g in gens(C)], false)
   l = preimage(ea,preimage(fa,d))
   return prod([C[i]^Int(l[i]) for i in 1:ngens(C)])
end

"""
    pol_elementary_divisors(x::MatElem)
    pol_elementary_divisors(x::MatrixGroupElem)

Return a list of pairs `(f_i,m_i)`, for irreducible polynomials `f_i` and
positive integers `m_i`, where the `f_i^m_i` are the elementary divisors of
`x`.
"""
function pol_elementary_divisors(A::MatElem{T}) where T
   a,_,c = _rational_canonical_form_setup(A)
   L = refine_for_jordan(a,c,A)[2]
   V = Vector{Tuple{typeof(L[1][2]),Int64}}(undef, length(L))
   for i in 1:length(L)
      V[i] = (L[i][2],L[i][3])
   end

# sorting the vector in order to have all equal polynomials in a consecutive bunch
# and block size in increasing order
   for i in 1:length(L)-1, j in i+1:length(L)
      if V[i][1]==V[j][1]
         V[i+1],V[j] = V[j],V[i+1]
         if V[i+1][2]<V[i][2] V[i],V[i+1] = V[i+1],V[i]  end
      end
   end

   return V
end

pol_elementary_divisors(x::MatrixGroupElem) = pol_elementary_divisors(x.elm)

"""
    generalized_jordan_block(f::T, n::Int) where T<:PolyElem

Return the Jordan block of dimension `n` corresponding to the polynomial `f`.
"""
function generalized_jordan_block(f::T, n::Int) where T<:PolyElem
   d = degree(f)
   JB = cat([companion_matrix(f) for i in 1:n]..., dims=(1,2))
   pos = 1
   for i in 1:n-1
      JB[pos:pos-1+degree(f), pos+degree(f):pos-1+2*degree(f)] = identity_matrix(base_ring(f),degree(f))
      pos += degree(f)
   end
   return JB
end

# TODO is there a way to accelerate the process? pol_elementary_divisors and generalized_jordan_block repeat parts of the same code.
"""
    generalized_jordan_form(A::MatElem{T}; with_pol::Bool=false) where T

Return (`J`,`Z`), where `Z^-1*J*Z = A` and `J` is a diagonal join of Jordan
blocks (corresponding to irreducible polynomials).
"""
function generalized_jordan_form(A::MatElem{T}; with_pol::Bool=false) where T
   V = pol_elementary_divisors(A)
   GJ = cat([generalized_jordan_block(v[1],v[2]) for v in V]..., dims=(1,2))
   a = rational_canonical_form(A)[2]
   gj = rational_canonical_form(GJ)[2]
   if with_pol return GJ, gj^-1*a, V
   else return GJ, gj^-1*a
   end
end


function is_conjugate(G::MatrixGroup, x::MatrixGroupElem, y::MatrixGroupElem)
   isdefined(G,:descr) || throw(ArgumentError("Group must be general or special linear group"))
   if G.descr==:GL || G.descr==:SL
      Jx,ax = jordan_normal_form(x.elm)
      Jy,ay = jordan_normal_form(y.elm)
      if Jx != Jy return false, nothing end
      z = inv(ax)*ay
      if G.descr==:GL return true, G(z) end
      ED = pol_elementary_divisors(x.elm)
      l = gcd([k[2] for k in ED])
      l = gcd(l, order(G.ring)-1)
      d = det(z)
      if isone(d^( div(order(G.ring)-1,l)))
         corr = _elem_given_det(x, d^-1)
         return true, G(corr*z)
      else return false, nothing
      end
   else
      return is_conjugate(G,x,y)
   end
end
