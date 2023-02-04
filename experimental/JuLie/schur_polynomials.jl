################################################################################
# Schur Polynomials
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Ulrich Thiel and OSCAR-ified by Claudia He Yun and Matthias Zach.
################################################################################

export schur_polynomial

"""
	schur_polynomial(λ::Partition{T}, n=sum(λ)::Int) where T<:Integer
	schur_polynomial(λ::Partition{T}, R::FmpzMPolyRing, n=sum(λ)::Int) where T<:Integer
	schur_polynomial(λ::Partition{T}, x::Array{fmpz_mpoly,1}) where T<:Integer

Returns the Schur polynomial ``s_λ(x₁,x₂,...,xₙ)`` in n variables, as a Multivariate Polynomial.

If neither `R` nor `x` are given, the Schur polynomial will be over `PolynomialRing(ZZ,["x1","x2",...,"xn"])`.

# Example
```julia-repl
julia> R,x = PolynomialRing(ZZ, ["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Integer Ring, fmpz_mpoly[a, b, c])
julia> schur_polynomial(Partition([2,1]),[x[1],x[2]])
a^2*b + a*b^2
julia> schur_polynomial(Partition([2,1]),R)
a^2*b + a^2*c + a*b^2 + 2*a*b*c + a*c^2 + b^2*c + b*c^2
julia> schur_polynomial(Partition([2]))
x1^2 + x1*x2 + x2^2
```

# Algorithm
We use two different Algorithms, depending on the size of the input.
The Combinatorial Algorithm is used for Partitions of small Integers, or if ``n ≥ 10``. In the other cases we use Cauchy's bialternant formula.

**Combinatorial Algorithm**
```math
s_λ:=∑_T x₁^{m₁}…xₙ^{mₙ}
```
where the sum is taken over all semistandard [tableaux](@ref JuLie.Tableau) ``T`` of shape ``λ``, and ``mᵢ`` gives the weight of ``i`` in ``T``.

**Cauchy's bialternant formula**
```math
s_λ(x₁,…,xₙ) =	∏_{1 ≤ i < j ≤ n} (x_i-x_j)^{-1} ⋅
\\begin{vmatrix}
x_1^{λ₁+n-1} & x_2^{λ_1+n-1} & … & x_n^{λ_1+n-1} \\\\
x_1^{λ_2+n-2} & x_2^{λ_2+n-2} & … & x_n^{λ_2+n-2} \\\\
⋮ & ⋮ & ⋱ & ⋮ \\\\
x_1^{λ_n} & x_2^{λ_n} & … & x_n^{λ_n}
\\end{vmatrix}
```
"""
function schur_polynomial(λ::Partition{T}, n=sum(λ)::Int) where T<:Integer
	n>=0 || throw(ArgumentError("n≥0 required"))
    if n==0 || n < length(λ)
        if isempty(λ)
            return 1
        else
            return 0
        end
    else
	    x = [string("x",string(i)) for i=1:n]
	    R,x = PolynomialRing(ZZ, x)
	    return schur_polynomial(R, λ, n)
    end
end


function schur_polynomial(R::FmpzMPolyRing, λ::Partition{T}, n=sum(λ)::Int) where T<:Integer
	n>=0 || throw(ArgumentError("n≥0 required"))
	if n > R.nvars
	n = R.nvars
	end
	x = gens(R)[1:n]
	if n==0 || n < length(λ)
	if isempty(λ)
		return 1
	else
		return 0
	end
	end

	if n>=10
	return schur_polynomial_combinat(R, λ, n)
	end
	#decide which Algorithm to use if n<10
	bo = sum(λ) <= [140,50,17,11,10,10,11,13,14][n]
	if bo
	return schur_polynomial_combinat(R, λ, n) #Combinatorial formula
	else
	return schur_polynomial_cbf(λ, x) #Cauchy's bialternant formula
	end
end


function schur_polynomial(λ::Partition{T}, x::Array{fmpz_mpoly,1}) where T<:Integer
	n = length(x)

	if n==0 || n < length(λ)
	if isempty(λ)
		return 1
	else
		return 0
	end
	end

	if n>=10
	R = x[1].parent
	return schur_polynomial_combinat(R, λ, n)
	end
	#decide which Algorithm to use if n<10
	bo = sum(λ) <= [140,50,17,11,10,10,11,13,14][n]
	if bo
	R = x[1].parent
	return schur_polynomial_combinat(R, λ, n) #Combinatorial formula
	else
	return schur_polynomial_cbf(λ, x) #Cauchy's bialternant formula
	end
end

#returning the schur polynomial in the first k generators of R using Cauchy's bialternant formula.
function schur_polynomial_cbf(λ::Partition{T}, x::Array{fmpz_mpoly,1}) where T<:Integer
	#if isempty(x) #this event is handled in the calling methods
	#	if sum(λ)==0
	#	return 1
	#	else
	#	return 0
	#	end
	#end

	n = length(x)
	R = x[1].parent # Multi-polynomialring
	S = R.base_ring # Integer Ring

	#if n < length(λ)
	#	return 0
	#end

	#=
	To calculate the determinant we use the Laplace expansion along the last row.
	Furthermore we use the fact, that each column consist of the same variable with decreasing powers.
	This allows us to factorize by the smallest power, thus our last row of the minors is always 1, this reduces the amount of polynomial multiplications we have to do.

	To avoid calculating the same minors multiple times, we calculate them from 1x1 up to nxn, always storing them in a Dictionary.
	Keep in mind that each minor consists of k columns and the top k rows.
	=#

	#initializing a few helpful Variables
	exponents = Int[getindex_safe(λ,i)+n-i for i=1:n] #the exponents from the Matrix read from top to bottom
	exp_incr = zeros(Int,n) #the increment with wich exponents increase
	for i = 1:n-1
	exp_incr[i] = exponents[i] - exponents[i+1]
	end
	exp_incr[n] = exponents[n]

	factor = one(R)
	d = R()

	#Initialize Dictionaries (calculating all possible combinations of k=1...n columns)
	sub_dets = [Dict{BitArray{1},fmpz_mpoly}() for i=1:n] #sub_dets[i] holds all the minors of size i, i.e. sub_dets[2][[false,true,true,false,false]] is the minor of the 2x2 matrix consisting of the first two rows intersected with the 2nd and 3rd columns.
	b = falses(n)
	ntrues = 0
	pointer = 1
	while true
	if b[pointer]
		while pointer<=n && b[pointer]
		b[pointer] = false
		ntrues -= 1
		pointer += 1
		end
		if pointer > n
		break
		end
	end
	b[pointer] = true
	ntrues += 1
	pointer = 1
	sub_dets[ntrues][copy(b)] = R()
	end

	#initialize sub_dets[1]
	for (columnview,) in sub_dets[1]
	sub_dets[1][columnview] = x[columnview][1]^(exp_incr[1])
	end

	#calculate sub_dets[2:n] using Laplace extension
	exp = zeros(Int,n)
	for i = 2:n
	for (columnview,) in sub_dets[i]
		d = R() #the alternating sum of minors
		s = ((i+n)%2==1 ? 1 : -1) #the alternating sign
		for j in (1:n)[columnview]
		columnview[j] = false
		d += s*sub_dets[i-1][columnview]
		columnview[j] = true
		s *= -1 #change sign
		end

		#multiply by the factorized term
		factor = MPolyBuildCtx(R)
		exp = zeros(Int,n)
		exp[columnview] .= exp_incr[i]
		push_term!(factor, one(S), exp)
		sub_dets[i][columnview] = mul!(d, d, finish(factor))
	end
	end
	sp = sub_dets[n][trues(n)]

	# divide by the product
	for i = 1:n-1
	for j = i+1:n
		sp = divexact(sp, x[i]-x[j])
	end
	end

	return sp
end

#returning the schur polynomial in the first k generators of R using the Combinatorial formula.
function schur_polynomial_combinat(R::FmpzMPolyRing, λ::Partition{T}, k=sum(λ)::Int) where T<:Integer
	if isempty(λ)
	return one(R)
	end

	S = R.base_ring
	sf = MPolyBuildCtx(R)

	#version of the function semistandard_tableaux(shape::Array{T,1}, max_val=sum(shape)::Integer)
	len = length(λ)
	Tab = [(fill(i,λ[i])) for i = 1:len]
	m = len
	n = λ[m]

	count = zeros(Int, R.nvars)
	valid = true
	while true
	count .= 0
	for i = 1:len
		for j = 1:λ[i]
		if Tab[i][j] <= k
			count[Tab[i][j]] += 1
		else
			valid = false
			break
		end
		end
	end
	if valid
		push_term!(sf, S(1), count)
	end
	#raise one element by 1
	while !(Tab[m][n] < k &&
		(n==λ[m] || Tab[m][n]<Tab[m][n+1]) &&
		(m==len || λ[m+1]<n || Tab[m][n]+1<Tab[m+1][n]))
		if n > 1
		n -= 1
		elseif m > 1
		m -= 1
		n = λ[m]
		else
		return finish(sf)
		end
	end

	Tab[m][n] += 1

	#minimize trailing elements
	if n < λ[m]
		i = m
		j = n + 1
	else
		i = m + 1
		j = 1
	end
	while (i<=len && j<=λ[i])
		if i == 1
		Tab[1][j] = Tab[1][j-1]
		elseif j == 1
		Tab[i][1] = Tab[i-1][1] + 1
		else
		Tab[i][j] = max(Tab[i][j-1], Tab[i-1][j] + 1)
		end
		if j < λ[i]
		j += 1
		else
		j = 1
		i += 1
		end
	end
	m = len
	n = λ[len]
	end #while true
end
