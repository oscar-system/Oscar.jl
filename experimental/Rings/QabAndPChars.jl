module QabModule

using Oscar
import Hecke: math_html
import Oscar: IJuliaMime
export QabField, QabAutomorphism, isconductor, root_of_unity, PChar, isroot_of_unity

###############################################################################
#
#   Implementation of QabField
#   Functions so that we can use Qab as coefficient field
#
###############################################################################

#
# Note that there are two possibilities construct a nth root of unity when n is even and n%4!=0
# either we can construct the field Q(z_n) or we take -z_(n/2) as a primitive n-th root
# to change between these two options, use PCharSaturateAll with allroots or allrootsNew (change this in
# the code)
#

#singular cannot use immutable structs
mutable struct QabField <: Nemo.Field # union of cyclotomic fields
end

function Base.show(io::IO, a::QabField)
  print(io, "Abelian closure of Q")
end

function Base.show(io::IO, ::IJuliaMime, a::QabField)
  io = IOContext(io, :compact => true)
  print(io, "\$")
  math_html(io, a)
  print(io, "\$")
end

function Hecke.math_html(io::IO, a::QabField)
  print(io, "\\text{Abelian closure of Q}")
end

mutable struct QabElem <: Nemo.FieldElem
  data::nf_elem  #actually in cyclotomic field
  c::Int #conductor of field
end

function Base.show(io::IO, a::QabElem)
  if get(io, :compact, false) == true
    math_html(io, a.data)
  else
    print(io, a.data, " in Q(z_$(a.c))")
  end
end

function Base.show(io::IO, ::IJuliaMime, a::QabElem)
  print(io, "\$")
  math_html(io, a)
  print(io, "\$")
end

function Hecke.math_html(io::IO, a::QabElem)
  if get(io, :compact, false) == true
    math_html(io, a.data)
  else
    print(io, a.data, " \\in Q(z_{$(a.c)})")
  end
end

function Oscar.singular_ring(F::QabField)
  return Singular.CoefficientRing(F)
end

function isconductor(n::Int)
  if isodd(n)
    return true
  end
  return n % 4 == 0
end

function coerce_up(K::AnticNumberField, n::Int, a::QabElem)
  d = div(n, a.c)
  @assert n % a.c == 0
  #z_n^(d) = z_a
  R = parent(parent(a.data).pol)
  return QabElem(evaluate(R(a.data), gen(K)^d), n)
end

function coerce_down(K::AnticNumberField, n::Int, a::QabElem)
  error("missing")
end

function make_compatible(a::QabElem, b::QabElem)
  if a.c == b.c
    return a,b
  end
  d = lcm(a.c, b.c)
  K, z = cyclotomic_field(d)
  return coerce_up(K, d, a), coerce_up(K, d, b)
end

import Base.+, Base.*, Base.-, Base.//, Base.==, Base.zero, Base.one, Base.^
import Nemo.mul!, Nemo.addeq!, Nemo.divexact, Nemo.iszero

Oscar.show_minus_one(::Type{QabElem}) = false
Oscar.needs_parentheses(a::QabElem) = Oscar.needs_parentheses(a.data)
Oscar.displayed_with_minus_in_front(a::QabElem) = Oscar.displayed_with_minus_in_front(a.data)

==(::QabField, ::QabField) = true

function ^(a::QabElem, n::Integer)
  return QabElem(a.data^n, a.c)
end

function ^(a::QabElem, n::fmpz)
  return a^Int(n)
end

function +(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data+b.data, a.c)
end

function addeq!(c::QabElem, a::QabElem)
	c, a=make_compatible(c, a)
	addeq!(c.data, a.data)
	return c
end

function -(a::QabElem)
  return QabElem(-a.data,a.c)
end

function neg!(a::QabElem)
	mul!(a.data,a.data,-1)
	return a
end

function *(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data*b.data, a.c)
end

*(a::Integer, b::QabElem) = QabElem(b.data*a, b.c)
*(a::fmpz, b::QabElem) = QabElem(b.data*a, b.c)

function mul!(c::QabElem, a::QabElem, b::QabElem)
	a, b = make_compatible(a, b)
	b, c = make_compatible(b, c)
  a, b = make_compatible(a, b)
	mul!(c.data, a.data, b.data)
	return c
end

function -(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data-b.data, a.c)
end

function //(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data//b.data, a.c)
end

# // with other name
function Base.div(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data//b.data, a.c)
end

function Oscar.divexact(a::QabElem, b::QabElem)
	a, b = make_compatible(a, b)
  return QabElem(divexact(a.data,b.data), a.c)
end

function Base.inv(a::QabElem)
	return(Base.one(Base.parent(a))//a)
end

function Oscar.isone(a::QabElem)
	return(isone(a.data))
end

function Base.iszero(a::QabElem)
	return(iszero(a.data))
end

import Base.==

function ==(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return a.data==b.data
end

function ==(a::QabElem, b::Integer)
	c = Base.parent(a)(b)
	a, c = make_compatible(a,c)
	return a==c
end

function (b::QabField)(a::Integer)
  return a*root_of_unity(b, 1)
end
function (b::QabField)(a::QabElem)
  return a
end
(b::QabField)() = b(0)
function (b::QabField)(a::Singular.n_unknown{QabElem}) 
  Singular.libSingular.julia(Singular.libSingular.cast_number_to_void(a.ptr))
end

function Base.copy(a::QabElem)
  return QabElem(a.data, a.c)
end

function Base.deepcopy(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data//b.data, a.c)
end

Base.parent(::QabElem) = QabField()
Base.one(::QabField) = QabField()(1)
Base.one(::QabElem) = QabField()(1)

Oscar.isnegative(::QabElem) = false

#perhaps the following is not necessary
#Oscar.promote_rule(::Type{QabElem}, ::Type{T}) where {T <: Integer} = QabElem

Oscar.promote_rule(::Type{QabElem}, ::Type{fmpz}) = QabElem

Oscar.promote_rule(::Type{QabElem}, ::Type{fmpq}) = QabElem

Oscar.promote_rule(::Type{QabElem}, ::Type{fmpq_poly}) = QabElem


###############################################################################
#
#   Functions for computing roots
#
###############################################################################

function root_of_unity(K::QabField, n::Int)
  #this function finds a primitive root of unity in our field, note this is not always e^(2*pi*i)/n
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(c)
  if c == n
    return QabElem(z, c)
  else
    return QabElem(-z, c)
  end
end

function root_of_unity2(K::QabField, n::Int)
  #this function returns the primitive root of unity e^(2*pi*i/n)
  c=n
  K, z = cyclotomic_field(c)
  return QabElem(z, c)
end

function Oscar.root(a::QabElem, n::Int)
 o = Oscar.order(a)
 l = o*n
 mu = root_of_unity2(QabField(), Int(l))
 return mu
end

function allroot(a::QabElem, n::Int)
	#all roots in a probably smaller field than with the function allrootNew
	#(using root_of_unity -> construct field Q(z_5) when needing a 10th root of unity,
	#root_of_unity constructs the field Q(z_10))
	o = Oscar.order(a)
  l = o*n
  mu = root_of_unity(QabField(), Int(l))
	A=QabElem[]
	if l==1 && mu^1==a
		A=[A;mu]
	end
  for k=1:(l-1)
		if (mu^k)^n==a
    	A=[A; mu^(k)]
		end
  end
  return A
end

function allrootNew(a::QabElem, n::Int)
	#compute all nth roots of a, where a has to be a root_of_unity
	o=Oscar.order(a)
	l=o*n
	A=QabElem[]
	mu=root_of_unity2(QabField(),Int(l))
	if l==1 && mu^1==a
		A=[A;mu]
	end
	for k=1:(l-1)
		if (mu^k)^n==a
			A=[A;mu^k]
		end
	end
	if size(A,1)==0
		error("no root found")
	end
	return A
end


###############################################################################
#
#   Galois automorphisms of Qab
#
#   The Galois automorphisms of the $n$-th cyclotomic field are the maps
#   defined by $\zeta_n \mapsto \zeta_n^k$, for $1 \leq k < n$,
#   with $\gcd(n, k) = 1$.
#   Thus we can define automorphisms $\sigma_k$ of Qab as follows.
#   For each prime power $q$, $\zeta_q$ is mapped to $\zeta_q^k$ if
#   $k$ and $q$ are coprime, and to $\zeta_q$ otherwise.
#
#   The action of such a map $\sigma_k$ on the $n$-th cyclotomic field can be
#   described by $\sigma_l$, with $l$ coprime to $n$:
#   Write $n = n_0 n_1$ where $\gcd(n_0, n_1) = 1$ and $n_1$ is maximal
#   with $\gcd(k, n_1) = 1$, and choose $a, b$ with $1 = a n_0 + b n_1$.
#   Then $l = k a n_0 + b n_1$ is coprime to $n$ and has the properties
#   $l \equiv 1 \pmod{n_0}$ and $l \equiv k \pmod{n_1}$.
#
###############################################################################

mutable struct QabAutomorphism
    exp::Int
end

function ^(val::QabElem, sigma::QabAutomorphism)
    k = sigma.exp
    n = val.c
    g = gcd(k, n)
    if g != 1
      # Replace `k` by an equivalent one that is coprime to `n`.
      n0 = 1
      n1 = n
      for (p, exp) in collect(Oscar.factor(g))
        while mod(n1, p) == 0
          n0 = n0*p
          n1 = div(n1, p )
        end
      end
      (gg, a, b) = gcdx(n0, n1)
      @assert gg == 1 "n0 and n1 should be coprime"
      k = k*a*n0 + b*n1
    end
    data = val.data  # nf_elem
    coeffs = Nemo.coefficients(data)
    res = zeros(eltype(coeffs), n)
    res[1] = coeffs[1]
    for i in 2:length(coeffs)
      res[mod((i-1)*k, n)+1] = coeffs[i]
    end
    F = parent(data) # cycl. field
    R = parent(F.pol)
    return QabElem(F(R(res)), n)
end


###############################################################################
#
#   Elements in quadratic subfields of cyclotomic fields
#
#
###############################################################################

function generators_galois_group_cyclotomic_field(n::Int)
    res = GAP.Globals.GeneratorsPrimeResidues(GAP.julia_to_gap(n))
    return [QabAutomorphism(k)
            for k in Vector{Int}(GAP.Globals.Flat(res.generators))]
end

"""
    square_root_in_cyclotomic_field(F::QabField, n::Int, N::Int)

Return an element `a` in the field `F` that is represented w.r.t. the `N`-th
cyclotomic field and has the property `a^2 == n`.

If the `N`-th cyclotomic field does not contain such an element
then `nothing` is returned.

If `n` is positive then `a` is the positive square root of `n`,
otherwise `a` is a positive multiple of the imaginary unit.
(Here we assume that the underlying primitive `N`-th root of unity
is identified with the complex number `exp(2*Pi*i/N)`,
where `i` is the imaginary unit.)
"""
function square_root_in_cyclotomic_field(F::QabField, n::Int, N::Int)
    z = root_of_unity(F, N)
    cf = 1
    sqf = 1
    for (p,e) in collect(factor(n))
      cf = cf * p^div(e, 2)
      if e % 2 != 0
        sqf = sqf * p
      end
    end
    nn = Int(sqf)  # nn is positive and squarefree
    @assert N % nn == 0

    n4 = nn % 4
    if n4 == 1
      if n < 0
        if N % 4 != 0
          return
        end
        N4 = div(N, 4)
        z4 = z^N4
        cf = cf * z4
      end
    elseif n4 == 2
      if N % 8 != 0
        return
      end
      N8 = div(N, 8)
      z8 = z^N8
      cf = cf * (z8 - z8^3)
      if n < 0
        cf = cf * z8^2
      end
    elseif n4 == 3
      if n > 0
        if N % 4 != 0
          return
        end
        N4 = div(N, 4)
        z4 = z^N4
        cf = cf * (-z4)
      end
    end

    # Compute the coefficients of the Atlas irrationality 2*b_nn+1,
    # w.r.t. the N-th cyclotomic field.
    # (The underlying formula is due to a theorem of Gauss.)
    cfs = zeros(fmpz, N)
    cfs[1] = 1
    q = div(N, nn)
    for k in 1:div(nn,2)
      pos = q * mod(k^2, nn) + 1
      cfs[pos] = cfs[pos] + 2
    end

    # Create the corresponding number field element.
    FF = parent(z.data)
    pol = FF.pol
    R = parent(pol)
    elm = mod(R(cfs), pol)

    return cf * Main.QabModule.QabElem(FF(elm), N)
end

"""
    quadratic_irrationality_info(a::QabModule.QabElem)

Return `(x, y, n)`, where `x`, `y` are of type `fmpq` and `n` is
a squarefree integer, such that `a == x + y sqrt(n)` holds.

(We assume that the underlying primitive `N`-th root of unity that
is used to define `a` is identified with the complex number `exp(2*Pi*i/N)`,
where `i` is the imaginary unit.)
"""
function quadratic_irrationality_info(a::QabModule.QabElem)
    n = a.c

    # Compute the Galois group generators of the `n`-th cyclotomic field.
    galgens = generators_galois_group_cyclotomic_field(n)

    # Start computing the orbit of `a` under the Galois group:
    # If the orbit length is larger than 2 then `a` is not in a
    # quadratic field, and we return nothing.
    cand = nothing
    for sigma in galgens
      img = a^sigma
      if img != a
        if isnothing(cand)
          cand = img
        elseif cand != img
          return
        end
      end
    end
    for sigma in galgens
      img = cand^sigma
      if img != a && img != cand
        return
      end
    end

    # We have a = x + y \sqrt{m} and cand = x - y \sqrt{m}.
    x = coeff(a.data + cand.data, 0) // 2
    root_multiple = a.data - x
    ysquarem = coeff(root_multiple^2, 0)  # fmpq
    num = numerator(ysquarem)
    den = denominator(ysquarem)
    den_y = sqrt(den)
    m = sign(num)
    for (p, e) in collect(factor(num))
      if e % 2 == 1
        m = m * p
      end
    end
    y = sqrt(ysquarem // m)

    # It remains to compute the sign of y.
    # (This relies on the choice of the primitive n-th root of unity
    # as the complex number exp(2*pi*i/n).)
    y_std_m = y * square_root_in_cyclotomic_field(parent(a), Int(m), n).data
    if y_std_m != root_multiple
      @assert y_std_m == - root_multiple
      y = -y
    end

    return (x, y, m)
end


###############################################################################
#
#   Partial character functions
#
###############################################################################

struct PChar
	#A has generators of the lattice in rows
  A::fmpz_mat
	#images of the generators are saved in b
  b::Array{FieldElem, 1}
	#Delta are the indices of the cellular variables of the associated ideal
	#(the partial character is a partial character on Z^Delta)
  D::Set{Int64}
end

function (Chi::PChar)(b::fmpz_mat)
  @assert Nemo.nrows(b)==1
  @assert Nemo.ncols(b) == Nemo.ncols(Chi.A)
  s = solve(Chi.A', b')
  return evaluate(FacElem(Dict([(Chi.b[i], s[i, 1]) for i=1:length(Chi.b)])))
end

function (Chi::PChar)(b::Array{Nemo.fmpz,})
	@assert size(b,2)==Nemo.ncols(Chi.A)
	B=Matrix(FlintZZ,1,Nemo.ncols(Chi.A),b)
	s = solve(Chi.A', B')
	return evaluate(FacElem(Dict([(Chi.b[i], s[i, 1]) for i=1:length(Chi.b)])))
end

function isroot_of_unity(a::QabElem)
  b = a^a.c
  return b.data == 1 || b.data == -1
end

function LatticeEqual(A::fmpz_mat,B::fmpz_mat)
  @assert Nemo.ncols(A)==Nemo.ncols(B)
  A=A'
  B=B'
  #use solve to check if lattices are equal
  testVector=matrix(FlintZZ,Nemo.nrows(A),1,zeros(Int64,Nemo.nrows(A),1))
  #test if A contained in B
  for k=1:Nemo.ncols(A)
		for j=1:Nemo.nrows(A)
			testVector[j,1]=A[j,k]
		end
		if cansolve(B,testVector)[1]==false
       	   return(false)
		end
  end

  for k=1:Nemo.ncols(A)
		for j=1:Nemo.nrows(A)
			testVector[j,1]=B[j,k]
		end
		if cansolve(A,testVector)[1]==false
      return(false)
		end
  end
  return(true)
end

function Oscar.order(a::QabElem)
  f = Nemo.factor(fmpz(2*a.c))
  o = 1
  for (p, e) = f.fac
    b = a^div(2*a.c, Int(p)^e)
    f = 0

    while !isone(b)
      b = b^p
      f += 1
    end
    o *= p^f
  end
  return o
end

function PCharEqual(P::PChar,Q::PChar)
  if LatticeEqual(P.A,Q.A)==false
		return false
  end

  #now test if the values taken on the generators of the lattices are equal
	for i=1:Nemo.nrows(P.A)
		TestVec=sub(P.A,i:i,1:Nemo.ncols(P.A))
		if P(TestVec)!=Q(TestVec)
			return false
		end
	end

	for i=1:Nemo.nrows(Q.A)
		TestVec=sub(Q.A,i:i,1:Nemo.ncols(P.A))
		if P(TestVec)!=Q(TestVec)
			return false
		end
	end
  return true
end

function Hecke.saturate(L::PChar)
  #this function doesn't work, we have to change root_of_unity to root_of_unity2 in the function root
  H = hnf(L.A')
  s = sub(H, 1:Nemo.ncols(H), 1:Nemo.ncols(H))
  i, d = pseudo_inv(s)  #is = d I_n
  #so, saturation is i' * H // d
  S = divexact(i'*L.A, d)
  l = QabElem[]
  for k=1:Nemo.nrows(s)
    c = i[1,k]
    for j=2:Nemo.ncols(s)
      c = gcd(c, i[j,k])
      if isone(c)
        break
      end
    end
    mu = evaluate(FacElem(Dict([(L.b[j], div(i[j, k], c)) for j=1:Nemo.ncols(s)])))
    mu = root(mu, Int(div(d, c)))
    push!(l,  mu)  # for all saturations, use all roots - a cartesian product
  end
  #new values are d-th root of l
  return PChar(S, l, L.D)
end

Oscar.elem_type(::QabField) = QabElem
Oscar.parent_type(::Type{QabElem}) = QabField
Oscar.parent_type(::QabElem) = QabField

function PCharSaturateAll(L::PChar)
	#computes all saturations of the partial character L
  Result=PChar[]

  #first handle case wher the domain of the partial character is the zero lattice
  #in this case return L
  ZeroTest=matrix(FlintZZ,1,Nemo.ncols(L.A),zeros(Int64,1,Nemo.ncols(L.A)))
  if LatticeEqual(L.A,ZeroTest)==true
		push!(Result,L)
		return Result
  end

  #now not trivial case
  H = hnf(L.A')
  s = sub(H, 1:Nemo.ncols(H), 1:Nemo.ncols(H))
  i, d = pseudo_inv(s)  #is = d I_n
  #so, saturation is i' * H // d
  S = divexact(i'*L.A, d)
  Re = QabElem[]

	B = Array[]
  for k=1:Nemo.nrows(s)
    c = i[1,k]
    for j=2:Nemo.ncols(s)
      c = gcd(c, i[j,k])
      if isone(c)
        break
      end
    end
    mu = evaluate(FacElem(Dict([(L.b[j], div(i[j, k], c)) for j=1:Nemo.ncols(s)])))
		#change between the two options allrootNew and allroot in order to use a possibly smaller or
		#bigger field
		#mu=allrootNew(mu, Int(div(d,c)))
		mu=allroot(mu, Int(div(d,c)))
    push!(B,  mu)
  end
  C=my_product(B)
  T=Array[]
  for a in C
		push!(T,collect(a))
  end

  for k=1:size(T,1)
		#check if PChar(S,T[k],L.D) puts on the right value on the lattice generators of L
		Pnew=PChar(S,T[k],L.D)
		flag=true	#flag if value on lattice generators is right
		for i=1:Nemo.nrows(L.A)
			if Pnew(sub(L.A,i:i,1:Nemo.ncols(L.A)))!=L.b[i]
				flag=false
				println("found wrong saturation (for information), we delete it")
			end
		end
		if flag==true
			push!(Result,PChar(S, T[k],L.D))
		end
  end
  return Result
end

function my_product(P::Array)
  T = ntuple(x->P[x], length(P))
  return Iterators.product(T...)
end

end
