##############################################################
#
# Conversion from Singular (weierstr.lib)
#
#############################################################

export weierstrass_preparation_theorem

export weierstrass_division_theorem


"""
    order_of_generality(f::MPolyElem, n::Int)

Compute for a given polynomial `f` an integer `d` such that `f` is `x_n`-general of order `d` 
with `x_n` the `n`-th variable of the polynomial ring of `f`.
"""
function order_of_generality(f::MPolyElem, n::Int)
	 
	vect = fill(parent(f)(0), nvars(parent(f))) 
	vect[n] = gen(parent(f),n) 
	f = evaluate(f, vect) 
	f != 0 || return false, -1
	return true, degree(leading_monomial(f, ordering=neglex(gens(parent(f)))), n)
end

"""
    jet(f::MPolyElem, degree_f::Int)

Given a polynomial `f` and an integer degree `degree_f` , this function returns the polynomial,
which contains only the monomials of `f`, which total degree is lower than or equal to `degree_f`.
"""
function jet(f::MPolyElem, degree_f::Int) 
  i = length(f)
  while i >= 1
    if total_degree(term(f,i)) > degree_f
      f = f - term(f,i)
    end
    i = i - 1
  end 
	return f
end


"""
    jet(f::MPolyElem, degree_f::Int, vecw::Vector{Int})

Given a polynomial `f`, an integer degree `degree_f` and an integer vector `vecw`, which has to have 
the same length as length(parent(f)), this function returns a polynomial, which only contains 
monomial, whose exponents taken as a vector multiplied (scalar product) with the given 
weight-vector `intw` are lower than or equal to `degree_f`.
"""
function jet(f::MPolyElem, degree_f::Int, vecw::Vector{Int}) 

	vect = collect(exponent_vectors(f))
	i = length(f)
	while i >= 1
		if transpose(vecw)*vect[i] > degree_f
			f = f-term(f, i)
		end
		i = i-1	
	end
	return f
end

"""
    invert_unit(u::MPolyElem, n::Int)
If `u` is an unit, returning boolean value for success and the inverse of `u`. 
Else return false and the zero polynomial. 
"""
function invert_unit(u::MPolyElem, n::Int)
	if u != 0 && total_degree(leading_monomial(u, ordering=negdegrevlex(gens(parent(u))))) == 0
		u0 = jet(u, 0)
		u = jet(1-divexact(u, u0), n)
		ui = u
		v = 1+u
		if u != 0
			i = div(n, total_degree(leading_monomial(u, ordering=negdegrevlex(gens(parent(u))))))
			while i > 1
				ui = jet(ui*u, n)
				v = v+ui
				i = i-1
			end
		end
		v = divexact(jet(v, n), u0)
		return true, v
	else
		return false, parent(u)(0)
	end
end

"""
    last_var_general(f::MPolyElem, n::Int)

Transform `f` into a polynomial, which is x_n general of finite order for x_n the `n`-th ring variable.
"""
var_general(f::MPolyElem, n::Int) = var_general(f, n, 0)

"""
    var_general(f::MPolyElem, n::Int, nesting_depth::Int)

Transform `f` into a polynomial, which is x_n general of finite order for x_n the `n`-th ring variable.
`nesting_depth` is the counter of recursive calls. First call with `nesting_depth` = 0.
"""
function var_general(f::MPolyElem, n::Int, nesting_depth::Int) 
	R = parent(f)
	j = nvars(R)
	_, d = order_of_generality(f, n) 
	if d > 0 
		return f	
  end 
	m = gens(R) 
	g = initial_form(f)
	i = 1
	while i <= j 
		if i != n
			m_hilf = m
			m_hilf[i] = R(0)
			if length(g) > length(evaluate(g, m_hilf))
				m[i] = gen(R,i)+rand(1-(nesting_depth)*10 : 1+(nesting_depth)*10)*gen(R,n)
				g = evaluate(f, m)
				break
			end
		end
		i = i+1
		end
		
		if nesting_depth <= 3
			return var_general(g, n, nesting_depth+1)
		end
		if nesting_depth == 4
			m=rand_coeff(m, n, 1, 1000, *)
		else
			m=rand_coeff(m, n, 2, nesting_depth*d, ^)
		end
		g = evaluate(f, m)
		return var_general(g, n, nesting_depth+1)		
end


"""
    initial_form(f::MPolyElem) 

Transform `f` into an initial form for the algorithm.
"""
function initial_form(f::MPolyElem)
	degree_f = total_degree(leading_monomial(f, ordering=negdegrevlex(gens(parent(f))))) 
	g = jet(f, degree_f) 
	return g
end

"""
    rand_coeff(m::Vector{fmpq_mpoly}, n::Int, mini::Int, maxi::Int, operator)

Given a vector of monomials `m`, an Integer `n`, integer bounds mini und maxi and an 
operator, this function computes the monomials for var_general, to compute a polynomial,
which is general for the `n`-th ring variable.
"""
function rand_coeff(m::Vector{fmpq_mpoly}, n::Int, mini::Int, maxi::Int, operator)
	R=parent(m[1])
	j=nvars(R)
	i = 1
	while i <= j 
		if i != n
			m[i] = gen(R,i)+operator(gen(R,n),rand(mini : maxi))
		end				
		i = i+1
	end
	return m
end 


"""
    weierstrass_division_theorem(g::MPolyElem, f::MPolyElem, degree::Int)

Given a multivariate polynomial `f`, which is general of order `d` in the last ring variable,
a polynomial `g` and an integer degree `degree`, performing the weierstrass division theorem of 
`f` up to degree `degree`, by dividing `f` by `g`. Returning a polynomial `u` and a remainder `r` 
up to degree `degree`, s.t. f=ug+r. Additional the number of iterations will be returned. 
"""
function weierstrass_division_theorem(g::MPolyElem, f::MPolyElem, degree::Int, n::Int)
	R = parent(f)
	r = R(0)
	vect = fill(0, nvars(R))
	vect[n] = 1 
	_, d = order_of_generality(f, n) 
	
	D = degree+d
	fhat = jet(f, d-1, vect) 
	ftilde = divexact((f-fhat), gen(R,n)^d)
	_, u = invert_unit(ftilde, D) 
	
	j = g 
	jhat = jet(j, d-1, vect) 
	jtilde = divexact(j-r, gen(R,n)^d)
	r = jhat
	h = jtilde
	i = 0
	while length(j) > 0

		j = jet(-fhat*u*jtilde, D) 
		jhat = jet(j, d-1, vect)
		jtilde = divexact(j-jhat, gen(R,n)^d)
		r = r+jhat
		h = h+jtilde
		i = i+1
	
	end
	return jet(u*h, degree), jet(r, degree), i

end

"""
    weierstrass_preparation_theorem(f::MPolyElem, degree::Int, n::Int)

Given a multivariate polynomial `f` and an integer degree `degree`, performing the 
weierstrass preparation theorem of `f` up to degree `degree` for the `n`-th ring variable.
Returning an unit `u` and a weierstrass polynomial `fhat`, s.t. `u``fhat`=`f`, 
while `fhat` is weierstrass. 
Additional the number of iterations used will be returned. 
"""
function weierstrass_preparation_theorem(f::MPolyElem, degree::Int, n::Int)
	
	if f == 0
		throw(DomainError(f, "polynomial is zero"))
	end
	R = parent(f)
	success, d = order_of_generality(f, n)
	if !success
		f = var_general(f, n, 0) 
		success, d = order_of_generality(f, n)
	end
	u, f2, iter_count = weierstrass_division_theorem(gen(R,n)^d, f, degree, n) 
	return u, gen(R,n)^d-f2, iter_count

end



