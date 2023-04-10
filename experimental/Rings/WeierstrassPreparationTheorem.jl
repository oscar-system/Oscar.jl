##############################################################
#
#Conversion from Singular (weierstr.lib)
#
#############################################################

"""
`order_of_generality(f::MPolyElem)`

compute for a given polynomial `f` an integer `d` such that `f` is `x_n`-general of order `d` with `x_n` the n-th variable of the polynomial ring of `f` 
"""
function order_of_generality(f::MPolyElem,n::Int)
	 
	vect=fill(parent(f)(0),nvars(parent(f))) 
	vect[n]=gens(parent(f))[n] 
	f=evaluate(f,vect) 
	f != 0 || return false,-1
	return true,degree(leading_monomial(f,ordering=neglex(gens(parent(f)))),n)
end

"""
`jet(f::MPolyElem,de::Int)`

Given a polynomial `f` and an integer degree `de` , this function returns the polynomial, which contains only the monomials of `f`, which total degree is lower than or equal to `de`
"""
function jet(f::MPolyElem,de::Int) 
	while f!=0
		if total_degree(leading_monomial(f,ordering=deglex(gens(parent(f)))))>de
			f=f-coeff(f,leading_monomial(f,ordering=deglex(gens(parent(f)))))*leading_monomial(f,ordering=deglex(gens(parent(f))))
		else 
			break
		end
	end
	return f
end


"""
`jet(f::MPolyElem,de::Int,vecw::Vector{Int})` 

Given a polynomial `f`, an integer degree `de` and an integer vector `vecw`, which has to have the same length as length(parent(f)), this function returns a polynomial, which only contains monomial, whose exponents taken as a vector multiplied (scalar product) with the given weight-vector `intw` are lower than or equal to `de`
"""
function jet(f::MPolyElem,de::Int,vecw::Vector{Int}) 

	vect=collect(exponent_vectors(f))
	i=length(f)
	while i>=1
		if transpose(vecw)*vect[i]>de
			f=f-term(f,i)
		end
		i=i-1	
	end
	return f
end

"""
`invert_unit(u::MPolyElem,n::Int)`
If `u` is an unit, returning boolean value for success and the inverse of `u`. Else return false and the zero polynomial. 
"""
function invert_unit(u::MPolyElem,n::Int)
	if u!=0 && total_degree(leading_monomial(u,ordering=negdegrevlex(gens(parent(u)))))==0
		u0=jet(u,0)
		u=jet(1-divexact(u,u0),n)
		ui=u
		v=1+u
		if u!=0
			i=div(n,total_degree(leading_monomial(u,ordering=negdegrevlex(gens(parent(u))))))
			while i > 1
				ui=jet(ui*u,n)
				v=v+ui
				i=i-1
			end
		end
		v=divexact(jet(v,n),u0)
		return true,v
	else
		return false, parent(u)(0)
	end
end

"""
`last_var_general(f::MPolyElem,n::Int)`

transform `f` into a polynomial, which is x_n general of finite order for x_n the `n`-th ring variable.
"""
function var_general(f::MPolyElem,n::Int)
	return var_general(f,n,0)
end

"""
`var_general(f::MPolyElem,n::Int,voice::Int)`

transform `f` into a polynomial, which is x_n general of finite order for x_n the `n`-th ring variable. `voice` is the counter of recursive calls (like nesting_depth). First call with `voice`=0
"""
function var_general(f::MPolyElem,n::Int,voice::Int) 
	R=parent(f)
	j=nvars(R)
	success,d=order_of_generality(f,n) 
	if d>0 
		return f	
	else 
		m=gens(R) 
		g=initial_form(f)
		i=1
		while i<=j 
			if i!=n
				m_hilf=m
				m_hilf[i]=R(0)
				if length(g) > length(evaluate(g,m_hilf))
					m[i]=gens(R)[i]+rand(1-(voice)*10:1+(voice)*10)*gens(R)[n]
					g=evaluate(f,m)
					break
				end
			end
			i=i+1
		end
		
		if voice <= 3
			return var_general(g,n,voice+1)
		end
		if voice == 4
			i=1
			while i<=j 
				if i!=n
					m[i]=gens(R)[i]+gens(R)[n]*rand(1:1000)
				end				
				i=i+1
			end 
			g=evaluate(f,m)
			return var_general(g,n,voice+1)
		else
			i=1
			while i<=j 
				if i!=n
					m[i]=gens(R)[i]+gens(R)[n]^rand(2:voice*d)
				end	
				i=i+1
			end
			g=evaluate(f,m)
			return var_general(g,voice+1)
		end
	end		
end


"""
`initial_form(f::MPolyElem)` 

transform `f` into an initial form for the algorithm
"""
function initial_form(f::MPolyElem)
	de=total_degree(leading_monomial(f,ordering=negdegrevlex(gens(parent(f))))) 
	g=jet(f,de) 
	return g
end

"""
`weierstrass_division_theorem(g::MPolyElem,f::MPolyElem,de::Int)`

Given a multivariate polynomial `f`, which is general of oder `d` in the last ring variable, a polynomial `g` and an integer degree `de`, performing the weierstrass division theorem of `f` up to degree `de`, by dividing `f` by `g`. Returning a polynomial `u` and a remainder `r` up to degree `de`, s.t. f=ug+r. Additional the number of iterations will be returned. 
"""
function weierstrass_division_theorem(g::MPolyElem,f::MPolyElem,de::Int,n::Int)
	R=parent(f)
	r=R(0)
	vect=fill(0,nvars(R))
	vect[n]=1 
	success,d=order_of_generality(f,n) 
	
	D=de+d
	fhat=jet(f,d-1,vect) 
	ftilde= divexact((f-fhat),gens(R)[n]^d)
	succ,u=invert_unit(ftilde,D) 
	
	j=g 
	jhat=jet(j,d-1,vect) 
	jtilde=divexact(j-r,gens(R)[n]^d)
	r=jhat
	h=jtilde
	i=0
	while length(j) >0

		j=jet(-fhat*u*jtilde,D) 
		jhat= jet(j,d-1,vect)
		jtilde=divexact(j-jhat,gens(R)[n]^d)
		r=r+jhat
		h=h+jtilde
		i=i+1
	
	end
	return jet(u*h,de),jet(r,de),i

end

"""
`weierstrass_preparation_theorem(f::MPolyElem,de::Int)`

Given a multivariate polynomial `f` and an integer degree `de`, performing the weierstrass preparation theorem of `f` up to degree `de`. Returning an unit `u` and a weierstrass polynomial `fhat`, s.t. `u``fhat`=`f`, while `fhat` is weierstrass. Additional the number of iterations used will be returned. 
"""
function weierstrass_preparation_theorem(f::MPolyElem,de::Int,n::Int)
	
	if f==0
		throw(DomainError(f, "polynomial is zero"))
	end
	R=parent(f)
	success,d=order_of_generality(f,n)
	if success!=true
		f=var_general(f,n,0) 
		success,d=order_of_generality(f,n)
	end
	u,f2,iter_count=weierstrass_division_theorem(gens(R)[n]^d,f,de,n) 
	return u,gens(R)[n]^d-f2,iter_count

end

