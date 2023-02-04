################################################################################
# Tableaux
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Ulrich Thiel and OSCAR-ified by Claudia He Yun and Matthias Zach.
################################################################################

export Tableau, shape, semistandard_tableaux, is_standard, is_semistandard, standard_tableaux, schensted, hook_length, hook_lengths, num_standard_tableaux, reading_word, weight, bump!

"""
	Tableau{T} <: AbstractArray{AbstractArray{T,1},1}

A **Young diagram** is a diagram of finitely many empty "boxes" arranged in left-justified rows, with the row lengths in non-increasing order. The box in row i and and column j has the **coordinates** (i,j). Listing the number of boxes in each row gives a partition λ of a non-negative integer n (the total number of boxes of the diagram). The diagram is then said to be of **shape** λ. Conversely, one can associate to any partition λ a Young diagram in the obvious way, so Young diagrams are just another way to look at partitions.

A **Young tableau** of shape λ is a filling of the boxes of the Young diagram of λ with elements from some set. After relabeling we can (and will) assume that we fill from a set of integers from 1 up to some number, which in applications is often equal to n. We encode a tableau as an array of arrays and we have implemented an own type ```Tableau{T}```	as subtype of ```AbstractArray{AbstractArray{T,1},1}``` to work with tableaux. As for partitions, you may increase performance by casting into smaller integer types, e.g.

For efficiency, we do not check whether the given array is really a tableau, i.e. whether the structure of the array defines a partition.

# Example
```julia-repl
julia> Tab=Tableau([[1,2,3],[4,5],[6]])
julia> Tab=Tableau(Array{Int8,1}[[2,1], [], [3,2,1]]) #Using 8 bit integers
```

# References
1. Wikipedia, [Young tableau](https://en.wikipedia.org/wiki/Young_tableau).
"""
struct Tableau{T} <: AbstractArray{AbstractArray{T,1},1}
	t::Array{Array{T,1},1}
end

function Base.show(io::IO, ::MIME"text/plain", Tab::Tableau)
	print(io, Tab.t)
end

function Base.size(Tab::Tableau)
	return size(Tab.t)
end

function Base.length(Tab::Tableau)
	return length(Tab.t)
end

function Base.getindex(Tab::Tableau, i::Int)
	return getindex(Tab.t,i)
end

function Base.getindex(Tab::Tableau, I::Vararg{Int, 2})
	return getindex(getindex(Tab.t,I[1]), I[2])
end




"""
	shape(Tab::Tableau{T})

Returns the shape of a tableau, i.e. the partition given by the lengths of the rows of the tableau.
"""
function shape(Tab::Tableau{T}) where T
	return Partition{T}([ length(Tab[i]) for i=1:length(Tab) ])
end


"""
	weight(Tab::Tableau)

The **weight** of a tableau is the number of times each number appears in the tableau. The return value is an array whose i-th element gives the number of times the integer i appears in the tableau.
"""
function weight(Tab::Tableau)
	if isempty(Tab)
		return Int[]
	end

	max = 0
	for i = 1:length(Tab)
		if max < Tab[i][end]
			max = Tab[i][end]
		end
	end

	w = zeros(Int,max)
	for rows in Tab
		for box in rows
			w[box] += 1
		end
	end
	return w
end


"""
	reading_word(Tab::Tableau)

The **reading word** of a tableau is the word obtained by concatenating the fillings of the rows, starting from the *bottom* row. The word is here returned as an array.

# Example
```
julia> reading_word(Tableau([ [1,2,3] , [4,5] , [6] ]))
6-element Array{Int64,1}:
 6
 4
 5
 1
 2
 3
```
"""
function reading_word(Tab::Tableau)
	w = zeros(Int,sum(shape(Tab)))
	k = 0
	for i = length(Tab):-1:1
		for j = 1:length(Tab[i])
			k += 1
			w[k] = Tab[i,j]
		end
	end
	return w
end


"""
	is_semistandard(Tab::Tableau)

A tableau is called **semistandard** if the entries weakly increase along each row and strictly increase down each column.
"""
function is_semistandard(Tab::Tableau)
	s = shape(Tab)
	if isempty(s)
		return true
	end

	#correct shape
	for i = 1:length(s)-1
		if s[i] < s[i+1]
			return false
		end
	end

	#increasing first row
	for j = 2:s[1]
		if Tab[1][j] < Tab[1][j-1]
			return false
		end
	end

	#increasing first column
	for i = 2:length(s)
		if Tab[i][1] <= Tab[i-1][1]
			return false
		end
	end

	#increasing rows and columns
	for i = 2:length(Tab)
		for j = 2:s[i]
			if Tab[i][j] < Tab[i][j-1] || Tab[i][j] <= Tab[i-1][j]
				return false
			end
		end
	end
	return true
end



"""
	semistandard_tableaux(shape::Partition{T}, max_val=sum(shape)::Integer) where T<:Integer

Returns a list of all semistandard tableaux of given shape and filling elements bounded by `max_val`. By default, `max_val` is equal to the sum of the shape partition (the number of boxes in the Young diagram). The list of tableaux is in lexicographic order from left to right and top to bottom.
"""
function semistandard_tableaux(shape::Partition{T}, max_val=sum(shape)::Integer) where T<:Integer
	SST = Array{Tableau{T},1}()
	len = length(shape)
	if max_val < len
		return SST
	elseif len==0
		push!(SST, Tableau(Array{T,1}[]))
		return SST
	end
	Tab = [Array{T}(fill(i,shape[i])) for i = 1:len]
	m = len
	n = shape[m]

	while true
		push!(SST,Tableau([copy(row) for row in Tab]))

		#raise one element by 1
		while !(Tab[m][n]<max_val &&
			(n==shape[m] || Tab[m][n]<Tab[m][n+1]) &&
			(m==len || shape[m+1]<n || Tab[m][n]+1<Tab[m+1][n]))
			if n > 1
				n -= 1
			elseif m > 1
				m -= 1
				n = shape[m]
			else
				return SST
			end
		end

		Tab[m][n] += 1

		#minimize trailing elements
		if n < shape[m]
			i = m
			j = n + 1
		else
			i = m + 1
			j = 1
		end
		while (i<=len && j<=shape[i])
			if i==1
				Tab[1][j] = Tab[1][j-1]
			elseif j==1
				Tab[i][1] = Tab[i-1][1] + 1
			else
				Tab[i][j] = max(Tab[i][j-1], Tab[i-1][j] + 1)
			end
			if j < shape[i]
				j += 1
			else
				j = 1
				i += 1
			end
		end
		m = len
		n = shape[len]
	end

end

"""
	semistandard_tableaux(shape::Partition{T}, max_val=sum(shape)::Integer) where T<:Integer

Shortcut for ```semistandard_tableaux(Partition(shape),max_val)```.
"""
function semistandard_tableaux(shape::Array{T,1}, max_val=sum(shape)::T) where T<:Integer
	return semistandard_tableaux(Partition(shape), max_val)
end

"""
	semistandard_tableaux(box_num::T, max_val=box_num::T) where T<:Integer

Returns a list of all semistandard tableaux consisting of `box_num` boxes and filling elements bounded by `max_val`.
"""
function semistandard_tableaux(box_num::T, max_val=box_num::T) where T<:Integer
	box_num>=0 || throw(ArgumentError("box_num ≥ 0 required"))
	SST = Array{Tableau{T},1}()
	if max_val<=0
		return SST
	end
	shapes = partitions(box_num)

	for s in shapes
		if max_val >= length(s)
			append!(SST, semistandard_tableaux(s.p,max_val))
		end
	end

	return SST
end


"""
	semistandard_tableaux(s::Array{T,1}, weight::Array{T,1}) where T<:Integer

Returns a list of all semistandard tableaux with shape s and given weight. This
requires that sum(s) = sum(weight).
"""
function semistandard_tableaux(s::Array{T,1}, weight::Array{T,1}) where T<:Integer
	n_max = sum(s)
	n_max==sum(weight) || throw(ArgumentError("sum(s) = sum(weight) required"))

	Tabs = Array{Tableau,1}()
	if isempty(s)
		push!(Tabs, Tableau(Array{Int,1}[]))
		return Tabs
	end
	ls = length(s)

	Tab = Tableau([ [0 for j = 1:s[i]] for i = 1:length(s)])
	sub_s = zeros(Integer, length(s))

	#tracker_row = zeros(Integer,n_max)

	function rec_sst!(n::Integer)

		#fill the remaining boxes if possible, else set them to 0
		if n == length(weight)
			for i = 1:ls
				for j = sub_s[i]+1:s[i]
					Tab[i][j] = n
					if i!=1 && Tab[i-1][j]==n
						for k = 1:i
							for l = sub_s[k]+1:s[k]
								Tab[i][j] = 0
							end
						end
						return
					end
				end
			end
			push!(Tabs,Tableau([copy(row) for row in Tab]))

			return

		#skip to next step if weight[n]==0
		elseif weight[n] == 0
			rec_sst!(n+1)
			return
		end

		#here starts the main part of the function
		tracker_row = zeros(Integer, weight[n])
		i = 1
		while sub_s[i] == s[i]
			i += 1
		end
		j = sub_s[i] + 1

		m=0
		while m >= 0
			if m == weight[n]		 #jump to next recursive step
				rec_sst!(n+1)
				Tab[tracker_row[m]][sub_s[tracker_row[m]]] = 0
				i = tracker_row[m] + 1
				if i <= ls
					j = sub_s[i] + 1
				end
				m -= 1
				sub_s[i-1] -= 1

			elseif i > ls
				if m == 0
					return
				else
					Tab[tracker_row[m]][sub_s[tracker_row[m]]] = 0
					i = tracker_row[m] + 1
					if i <= ls
						j = sub_s[i] + 1
					end
					m -= 1
					sub_s[i-1] -= 1
				end

			elseif j<=s[i] && (i==1 || (j<=sub_s[i-1] && n>Tab[i-1][j]))	#add an entry
				m += 1
				Tab[i][j] = n
				sub_s[i] += 1
				tracker_row[m] = i
				j += 1

			else #move pointerhead
				i += 1
				if i <= ls
					j = sub_s[i] + 1
				end
			end
		end	#while

	end	#rec_sst!()

	rec_sst!(1)
	return Tabs
end

function semistandard_tableaux(s::Partition{T}, weight::Partition{T}) where T<:Integer
	return semistandard_tableaux(Array{T,1}(s),Array{T,1}(weight))
end



"""
	is_standard(Tab::Tableau)

A tableau is called **standard** if it is semistandard and the entries are in bijection with 1,…,n, where n is the number of boxes.
"""
function is_standard(Tab::Tableau)
	s = shape(Tab)
	if isempty(s)
		return true
	end

	#correct shape
	for i = 1:length(s)-1
	if s[i] < s[i+1]
		return false
	end
	end

	#contains all numbers from 1 to n
	n = sum(s)
	numbs = falses(n)
	for i = 1:length(s)
		for j = 1:s[i]
			if Tab[i][j]>n
				return false
			end
			numbs[Tab[i][j]] = true
		end
	end
	if false in numbs
		return false
	end

	#increasing first row
	for j = 2:s[1]
		if Tab[1][j] <= Tab[1][j-1]
			return false
		end
	end

	#increasing first column
	for i = 2:length(s)
		if Tab[i][1] <= Tab[i-1][1]
			return false
		end
	end

	#increasing rows and columns
	for i = 2:length(s)
		for j = 2:s[i]
			if Tab[i][j] <= Tab[i][j-1] || Tab[i][j] <= Tab[i-1][j]
				return false
			end
		end
	end
	return true
end


"""
	standard_tableaux(s::Partition)
	standard_tableaux(s::Array{Integer,1})

Returns a list of all standard tableaux of a given shape.
"""
function standard_tableaux(s::Partition)
	Tabs = Array{Tableau,1}()
	if isempty(s)
		push!(Tabs, Tableau(Array{Int,1}[]))
		return Tabs
	end
	n_max = sum(s)
	ls = length(s)

	Tab = Tableau([ [0 for j = 1:s[i]] for i = 1:length(s)])
	sub_s = [0 for i=1:length(s)]
	Tab[1][1] = 1
	sub_s[1] = 1
	tracker_row = [0 for i=1:n_max]
	tracker_row[1] = 1

	n = 1
	i = 1
	j = 2

	while n > 0
		if n == n_max || i > ls
			if n == n_max
				push!(Tabs,Tableau([copy(row) for row in Tab]))
			end
			Tab[tracker_row[n]][sub_s[tracker_row[n]]] = 0
			i = tracker_row[n] + 1
			if i <= ls
				j = sub_s[i] + 1
			end
			n -= 1
			sub_s[i-1] -= 1
		elseif j<=s[i] && (i==1 || j<=sub_s[i-1])
			n += 1
			Tab[i][j] = n
			sub_s[i] += 1
			tracker_row[n] = i
			i = 1
			j = sub_s[1] + 1
		else
			i += 1
			if i <= ls
				j = sub_s[i] + 1
			end
		end
	end

	return Tabs
end


function standard_tableaux(s::Array{T,1}) where T<:Integer
	return standard_tableaux(Partition(s))
end


"""
	standard_tableaux(n::Integer)

Returns a list of all standard tableaux with n boxes.
"""
function standard_tableaux(n::Integer)
	n>=0 || throw(ArgumentError("n ≥ 0 required"))
	ST = Array{Tableau,1}()
	for s in partitions(n)
		append!(ST, standard_tableaux(s))
	end
	return ST
end



"""
	hook_length(lambda::Partition, i::Integer, j::Integer)

Consider the Young diagram of a partition λ. The **hook length** of a box, is the number of boxes to the right in the same row + the number of boxes below in the same column + 1. The function returns the hook length of the box with coordinates (i,j). The functions assumes that the box exists.
"""
function hook_length(lambda::Partition, i::Integer, j::Integer)
	h = lambda[i] - j + 1
	k = i + 1
	while k<=length(lambda) && lambda[k]>=j
		k += 1
		h += 1
	end
	return h
end

"""
	hook_length(Tab::Tableau, i::Integer, j::Integer)

Shortcut for ```hook_length(shape(Tab),i,j)```.
"""
function hook_length(Tab::Tableau, i::Integer, j::Integer)
	return hook_length(shape(Tab),i,j)
end


"""
	hook_lengths(lambda::Partition)

Returns the tableau of shape λ in which the entry at position (i,j) is equal to the hook length of the corresponding box.
"""
function hook_lengths(lambda::Partition)
	if isempty(lambda)
		return Tableau(Array{Int,1}[])
	end
	Tab = [ [hook_length(lambda,i,j) for j in 1:lambda[i]] for i in 1:length(lambda) ]
	return Tableau(Tab)
end




@doc raw"""
	num_standard_tableaux(lambda::Partition)

Returns the number ``f^\lambda`` of standard tableaux of shape ``λ`` using the hook length formula
```math
f^{\lambda} = \frac{n!}{\prod_{i,j} h_\lambda(i,j)} \;,
```
where the product is taken over all boxes in the Young diagram of ``\lambda`` and ``h_\lambda`` denotes the hook length of the box (i,j).

# References
1. Wikipedia, [Hook length formula](https://en.wikipedia.org/wiki/Hook_length_formula).
"""
function num_standard_tableaux(lambda::Partition)
	n = sum(lambda)
	h=factorial(ZZ(n))
	for i = 1:length(lambda)
		for j = 1:lambda[i]
			h = div(h,ZZ(hook_length(lambda,i,j)))
		end
	end
	return h
end


"""
	schensted(sigma::Array{Integer,1})
	schensted(sigma::Perm{T})

The Robinson–Schensted correspondence is a bijection between permutations and pairs of standard Young tableaux of the same shape. For a permutation sigma (given as an array), this function performs the Schnested algorithm and returns the corresponding pair of standard tableaux (the insertion and recording tableaux).

# Example
```julia-repl
julia> P,Q = schensted([3,1,6,2,5,4]);
julia> P
[[1, 2, 4], [3, 5], [6]]

julia> Q
[[1, 3, 5], [2, 4], [6]]
```

# References
1. Wikipedia, [Robinson–Schensted correspondence](https://en.wikipedia.org/wiki/Robinson–Schensted_correspondence)
"""
function schensted(sigma::Array{T,1}) where T<:Integer
	if isempty(sigma)
		return Tableau(Array{T,1}[]),Tableau(Array{T,1}[])
	end
	P = Tableau{T}([[sigma[1]]])
	Q = Tableau{T}([[1]])
	for i = 2:length(sigma)
		bump!(P, sigma[i], Q, i)
	end
	return P,Q
end

function schensted(sigma::Perm{T}) where T<:Integer
	return schensted(sigma.d)
end


"""
	bump!(Tab::Tableau, x::Int)

Inserts the integer x into the tableau Tab according to the bumping algorithm by applying the Schensted insertion.

# References
1. Wolfram MathWorld, [Bumping Algorithm](https://mathworld.wolfram.com/BumpingAlgorithm.html)
"""
function bump!(Tab::Tableau, x::Integer)
	if isempty(Tab)
		push!(Tab.t,[x])
		return Tab
	end

	i = 1
	while i <= length(Tab)
		if Tab[i, length(Tab[i])] <= x
			push!(Tab[i], x)
			return Tab
		end
		j = 1
		while j <= length(Tab[i])
			if Tab[i,j] > x
				temp = x
				x = Tab[i,j]
				Tab[i][j] = temp
				i += 1
				break
			end
			j += 1
		end
	end
	push!(Tab.t,[x])
	return Tab
end

"""
	bump!(Tab::Tableau, x::Integer, Q::Tableau, y::Integer)

Inserts x into Tab according to the bumping algorithm by applying the Schensted insertion. Traces the change with Q by inserting y at the same Position in Q as x in Tab.
"""
function bump!(Tab::Tableau, x::Integer, Q::Tableau, y::Integer)
	if isempty(Tab)
		push!(Tab.t,[x])
		push!(Q.t,[x])
		return Tab,Q
	end
	i = 1
	while i <= length(Tab)
		if Tab[i,length(Tab[i])] <= x
			push!(Tab[i], x)
			push!(Q[i], y)
			return Tab
		end
		j=1
		while j <= length(Tab[i])
			if Tab[i,j] > x
				temp = x
				x = Tab[i,j]
				Tab[i][j] = temp
				i += 1
				break
			end
		j += 1
		end
	end
	push!(Tab.t, [x])
	push!(Q.t, [y])

	return Tab,Q
end
