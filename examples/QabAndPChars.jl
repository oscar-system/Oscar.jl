module QabModule

using Oscar
import Hecke: math_html
import Oscar: IJuliaMime
export QabField, isconductor, root_of_unity, PChar, isroot_of_unity

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

function Oscar.isone(a::QabElem)
  return isone(a.data)
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
