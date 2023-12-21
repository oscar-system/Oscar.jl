using Oscar

global infoLevel = 0

###############################################################
# Implementation of different variants of the Groebner Walk.
# The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
###############################################################
@doc raw"""
		groebnerwalk(
			I::MPolyIdeal; 
			startOrder::MonomialOrdering = default_ordering(base_ring(I)),
			targetOrder::Symbol = lex(base_ring(I)),
			walktype::Symbol = :standard,
			perturbationDegree::Int = 2,
		)
Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using the Groebner Walk.
The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
One can choose a strategy of:
Standard Walk (:standard) computes the Walk like it´s presented in Cox, Little & O´Shea (2005).
Generic Walk (:generic) computes the Walk like it´s presented in Fukuda, Jensen, Lauritzen & Thomas (2005).
Perturbed Walk (:perturbed, with p = degree of the perturbation) computes the Walk like it´s presented in Amrhein, Gloor & Küchlin (1997).
Tran´s Walk (:tran) computes the Walk like it´s presented in Tran (2000).
Fractal Walk (:fractalcombined) computes the Walk like it´s presented in Amrhein & Gloor (1998) with multiple extensions. The target monomial order has to be lex. This version uses the Buchberger Algorithm to skip weight vectors with entries bigger than Int32.
Fractal Walk (:fractal) computes the Walk like it´s presented in Amrhein & Gloor (1998). Pertubes only the target vector.
Set infoLevel to trace the walk:
	-'infoLevel=0': no detailed printout,
	-'infoLevel=1': adds intermediate weight vectors,
	-'infoLevel=2': adds information about the Groebner basis.

#Arguments
- `G::Oscar.IdealGens`: ideal one wants to compute a Groebner basis for.
- `startOrder::MonomialOrdering=:degrevlex`: monomial order to begin the conversion.
- `targetOrder::MonomialOrdering=:lex`: monomial order one wants to compute a Groebner basis for.
- `walktype::Symbol=standard`: strategy of the Groebner Walk. One can choose a strategy of:
	- `standard`: Standard Walk,
	- `perturbed`: Perturbed Walk,
	- `tran`: Tran´s Walk,
	- `generic`: Generic Walk,
	- `fractal`: standard-version of the Fractal Walk,
	- `fractalcombined`: combined Version of the Fractal Walk. Target monomial order needs to be lex,
- `perturbationDegree::Int=2`: perturbationdegree for the perturbed Walk.


# Examples

```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"], ordering = :degrevlex);

julia> I = ideal([y^4+ x^3-x^2+x,x^4]);

julia> groebnerwalk(I, degrevlex(R), lex(R), :standard)
standard_walk results
Crossed Cones in: 
[4, 3]
[4, 1]
[12, 1]
[1, 0]
Cones crossed: 4
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])

julia> groebnerwalk(I, degrevlex(R), lex(R), :perturbed, 2)
perturbed_walk results
Crossed Cones in: 
[4, 3]
[4, 1]
[5, 1]
[12, 1]
[1, 0]
Cones crossed: 5
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])
```
"""
function groebnerwalk(
	I::MPolyIdeal, startOrder::MonomialOrdering = default_ordering(base_ring(I)),
	targetOrder::MonomialOrdering = lex(base_ring(I)),
	walktype::Symbol = :standard,
	perturbationDegree::Int = 2)
	S = canonical_matrix(startOrder)
	T = canonical_matrix(targetOrder)
	return groebnerwalk(
		I,
		S, T,
		walktype,
		perturbationDegree,
	)
end

@doc raw"""
		groebnerwalk(
			G::Oscar.IdealGens,
			S::Union{Matrix{T}, MatElem{T}},
			T::Union{Matrix{T}, MatElem{T}},
			walktype::Symbol = :standard,
			p::Int = 2
		) where T
Computes a reduced Groebner basis w.r.t. the monomial order T by converting the reduced Groebner basis G w.r.t. the monomial order S using the Groebner Walk.
One can choose a strategy of:
Standard Walk (:standard) computes the Walk like it´s presented in Cox et al. (2005).
Generic Walk (:generic) computes the Walk like it´s presented in Fukuda et al. (2005).
Perturbed Walk (:perturbed, with p = degree of the perturbation) computes the Walk like it´s presented in Amrhein et al. (1997).
Tran´s Walk (:tran) computes the Walk like it´s presented in Tran (2000).
Fractal Walk (:fractal) computes the Walk like it´s presented in Amrhein & Gloor (1998). Pertubes only the target vector.
Fractal Walk (:fractalcombined) computes the Walk like it´s presented in Amrhein & Gloor (1998) with multiple extensions. The target monomial order has to be lex. This version uses the Buchberger Algorithm to skip weightvectors with entries bigger than Int32.
Set infoLevel to trace the walk:
	-'infoLevel=0': no detailed printout,
	-'infoLevel=1': adds intermediate weight vectors,
	-'infoLevel=2': adds information about the Groebner basis.

#Arguments
-`G::Oscar.IdealGens`: Groebner basis to convert to the Groebner basis w.r.t. the target order T. G needs to be a Groebner basis w.r.t. the start order S.
-`S::Union{Matrix{N}, MatElem{N}}`: start monomial order w.r.t. the Groebner basis G. Note that S has to be a nxn-matrix with rank(S)=n and its first row needs to have positive entries.
-`T::Union{Matrix{N}, MatElem{N}}`: target monomial order one wants to compute a Groebner basis for. Note that T has to be a nxn-matrix with rank(T)=n and its first row needs to have positive entries.
-`walktype::Symbol=standard`: strategy of the Groebner Walk. One can choose a strategy of:
	- `standard`: Standard Walk (default),
	- `perturbed`: Perturbed Walk,
	- `tran`: Tran´s Walk,
	- `generic`: Generic Walk,
	- `fractal`: standard-version of the Fractal Walk,
	- `fractalcombined`: combined version of the Fractal Walk. The target monomial order needs to be lex,
-`p::Int=2`: perturbationdegree for the perturbed Walk.

# Examples

```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"], ordering = :degrevlex);

julia> I = ideal([y^4+ x^3-x^2+x,x^4]);

julia> groebnerwalk(I, [1 1; 0 -1], [1 0; 0 1], :standard)
standard_walk results
Crossed Cones in: 
[4, 3]
[4, 1]
[12, 1]
[1, 0]
Cones crossed: 4
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])

julia> groebnerwalk(I, [1 1; 0 -1], [1 0; 0 1], :perturbed, 2)
perturbed_walk results
Crossed Cones in: 
[4, 3]
[4, 1]
[5, 1]
[12, 1]
[1, 0]
Cones crossed: 5
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])
```
"""
function groebnerwalk(
	G::Union{Oscar.IdealGens, MPolyIdeal},
	S::Union{Matrix{N}, MatElem{N}},
	T::Union{Matrix{N}, MatElem{N}},
	walktype::Symbol = :standard,
	p::Int = 2,
) where N

	if walktype == :standard
		walk = (x) -> standard_walk(x, S, T)
	elseif walktype == :generic
		walk = (x) -> generic_walk(x, S, T)
	elseif walktype == :perturbed
		walk = (x) -> perturbed_walk(x, S, T, p)
	elseif walktype == :fractal
		walk = (x) -> fractal_walk(x, S, T)
	elseif walktype == :fractal_start_order
		walk = (x) -> fractal_walk_start_order(x, S, T)
		#    elseif walktype == :fractal_lex
		#       walk = (x) -> fractal_walk_lex(x, S, T)
		#   elseif walktype == :fractal_look_ahead
		#       walk = (x) -> fractal_walk_look_ahead(x, S, T)
	elseif walktype == :tran
		walk = (x) -> tran_walk(x, S, T)
	elseif walktype == :fractal_combined
		walk = (x) -> fractal_walk_combined(x, S, T)
	end
	delete_step_counter()
	S = Matrix{Int}(S)
	T = Matrix{Int}(T)
	# Make sure G is a fully reduced Groebner Basis
	Gb = groebner_basis(ideal(gens(G)), ordering = matrix_ordering(base_ring(G), S), complete_reduction = true)
	Gb = walk(Gb)

	if infoLevel >= 1
		println("Cones crossed: ", delete_step_counter())
	end
	return Oscar.IdealGens(gens(Gb), matrix_ordering(base_ring(G), T), isGB = true)
end

###########################################
# Counter for the steps
###########################################
counter = 0
function delete_step_counter()
	global counter
	temp = counter
	counter = 0
	return temp
end
function getstep_counter()
	global counter
	return counter
end
function raise_step_counter()
	global counter = getstep_counter() + 1
end
###############################################################
# Implementation of the standard walk.
###############################################################

function standard_walk(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
)
	if infoLevel >= 1
		println("standard_walk results")
		println("Crossed Cones in: ")
	end

	Gb = standard_walk(G, S, T, S[1, :], T[1, :])

	return Gb
end

function standard_walk(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
	currweight::Vector{Int},
	tarweight::Vector{Int})
	while true
		G = standard_step(G, currweight, T)
		if currweight == tarweight
			return G
		else
			currweight = next_weight(G, currweight, tarweight)
		end
		if infoLevel >= 1
			raise_step_counter()
			println(currweight)
			if infoLevel == 2
				println(G)
			end
		end
	end
end

###############################################################
# The standard step is used for the strategies standard and perturbed.
###############################################################

function standard_step(G::Oscar.IdealGens, w::Vector{Int}, T::Matrix{Int})
	R = base_ring(G)
	ordAlt = G.ord
	Gw = 0

	#check if no entry of w is bigger than Int32. If it´s bigger multiply it by 0.1 and round.
	if !checkInt32(w)
		Gw = ideal(initials(G, w))
		w, b = truncw(G, w, gens(Gw))
		!b && throw(error("Some entries of the intermediate weight-vector $w are bigger than int32"))
	else
		Gw = ideal(initials(G, w))
	end
	ordNew = create_order(R, w, T)
	H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)
	#H = liftGW2(G, ordAlt, Gw, H, ordNew)
	H = lift(G, ordAlt, H, ordNew)
	return interreduce_walk(H)
end

###############################################################
# Generic-version of the Groebner Walk.
###############################################################

function generic_walk(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
)
	Lm = [leading_term(g, ordering = matrix_ordering(base_ring(G), S)) for g in G]
	v = next_gamma(G, Lm, [0], S, T)
	ordNew = matrix_ordering(base_ring(G), T)
	if infoLevel >= 1
		println("generic_walk results")
		println("Facets crossed for: ")
	end
	while !isempty(v)
		G, Lm = generic_step(G, Lm, v, ordNew)
		raise_step_counter()

		if infoLevel >= 1
			println(v)
			if infoLevel == 2
				println(G)
			end
		end
		G = Oscar.IdealGens(G, ordNew)
		v = next_gamma(G, Lm, v, S, T)
	end
	return G
end

function generic_step(
	G::Oscar.IdealGens,
	Lm::Vector{T},
	v::Vector{Int},
	ord::MonomialOrdering,
) where {T <: MPolyRingElem}
	facet_Generators = facet_initials(G, Lm, v)
	H = groebner_basis(ideal(facet_Generators), ordering = ord, complete_reduction = true, algorithm = :buchberger)
	H, Lm = lift_generic(gens(G), Lm, gens(H), ord)
	G = interreduce(H, Lm, ord)
	return G, Lm
end

###############################################################
# Perturbed-version of the Groebner Walk.
###############################################################

function perturbed_walk(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
	p::Int,
)
	if infoLevel >= 1
		println("perturbed_walk results")
		println("Crossed Cones in: ")
	end

	currweight = perturbed_vector(G, S, p)

	while true
		tarweight = perturbed_vector(G, T, p)
		Tn = add_weight_vector(tarweight, T)
		G = standard_walk(G, S, Tn, currweight, tarweight)
		if same_cone(G, T)
			return G
		else
			p = p - 1
			currweight = tarweight
			S = Tn
		end
	end
end

###############################################################
# The Fractal Walk
###############################################################

##########################################
# global weightvectors
##########################################
pTargetWeights = []
pStartWeights = []
firstStepMode = false

###############################################################
# Combined version of the extensions of the Fractal Walk.
# This version
# - checks if the starting weight vector represents the monomial order and pertubes it if necessary.
# - analyses the Groebner basis Gw of the initialforms and uses the Buchberger-algorithm if the generators of Gw are binomial or less.
# - skips a step in top level in the last step.
# - checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################
function fractal_walk_combined(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
)
	if infoLevel >= 1
		println("fractal_walk_combined results")
		println("Crossed Cones in: ")
	end

	global pTargetWeights =
		[perturbed_vector(G, T, i) for i ∈ 1:nvars(base_ring(G))]
	return fractal_walk_combined(G, S, T, S[1, :], pTargetWeights, 1)
end

function fractal_walk_combined(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
	currweight::Vector{Int},
	pTargetWeights::Vector{Vector{Int}},
	p::Int,
)
	R = base_ring(G)
	G.isGB = true
	w = currweight
	ordAlt = G.ord

	# Handling the weight of the start order.
	if (p == 1)
		if !ismonomial(initials(G, w))
			global pStartWeights = [perturbed_vector(G, S, i) for i ∈ 1:nvars(R)]
			global firstStepMode = true
		end
	end
	if firstStepMode
		w = pStartWeights[p]
	else
		w = currweight
	end

	# main loop
	while true
		t = next_weight_fractal_walk(G, w, pTargetWeights[p])

		# Handling the final step in the current depth.
		# next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
		# -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.
		if t == 1 && p != 1
			if same_cone(G, T)
				if infoLevel >= 1
					println("depth $p: in cone ", currweight, ".")
				end

				# Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
				if !inCone(G, T, pTargetWeights, p)
					global pTargetWeights =
						[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
					if infoLevel >= 1
						println("depth $p: not in cone ", pTargetWeights[p], ".")
					end
				end
				return G
			end
		elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
			if inCone(G, T, pTargetWeights, p)
				if infoLevel >= 1
					println("depth $p: in cone ", pTargetWeights[p], ".")
				end
				return G
			end
			global pTargetWeights =
				[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
			if infoLevel >= 1
				println("depth $p: not in cone ", pTargetWeights[p], ".")
			end
			continue
		end

		# skip a step for target monomial order lex.
		if t == 1 && p == 1
			if infoLevel >= 1
				println("depth $p: recursive call in ", pTargetWeights[p])
			end
			return fractal_walk_combined(
				G,
				S,
				T,
				w,
				pTargetWeights,
				p + 1,
			)
		else
			w = w + t * (pTargetWeights[p] - w)
			w = convert_bounding_vector(w)
			Gw = ideal(initials(G, w))

			# handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
			if !checkInt32(w)
				println("depth $p: Some entries of $w are bigger than Int32. Trying to find another weight,")
				w, b = truncw(G, w, gens(Gw))
				if !b
					println("depth $p: Doing a direct conversion to the target monomial ordering.")
					ordNew = matrix_ordering(R, T)
					w = T[1, :]
					G = groebner_basis(ideal(G), ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

					if !inCone(G, T, pTargetWeights, p)
						global pTargetWeights =
							[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
						if infoLevel >= 1
							println("depth $p: not in cone ", pTargetWeights[p], ".")
						end
					end
					return G
				end
			end
			ordNew = create_order(R, w, T)
			# converting the Groebner basis
			if (p == nvars(R) || isbinomial(gens(Gw)))
				H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

				if infoLevel >= 1
					println("depth $p: conversion in ", w, ".")
				end
				raise_step_counter()
			else
				if infoLevel >= 1
					println("depth $p: recursive call in $w.")
				end
				H = fractal_walk_combined(
					Oscar.IdealGens(R, gens(Gw), ordAlt),
					S,
					T,
					deepcopy(currweight),
					pTargetWeights,
					p + 1,
				)
				global firstStepMode = false
			end
		end
		#H = liftGW2(G, ordAlt, Gw,H, ordNew)
		H = lift_fractal_walk(G, H, ordNew)
		G = interreduce_walk(H)
		ordAlt = ordNew
		currweight = w
	end
end

###############################################################
# Plain version of the Fractal Walk.
# This version checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################

function fractal_walk(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
)
	if infoLevel >= 1
		println("FractalWalk_standard results")
		println("Crossed Cones in: ")
	end

	global pTargetWeights =
		[perturbed_vector(G, T, i) for i ∈ 1:nvars(base_ring(G))]
	return fractal_recursiv(G, S, T, S[1, :], pTargetWeights, 1)
end

function fractal_recursiv(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
	currweight::Vector{Int},
	pTargetWeights::Vector{Vector{Int}},
	p::Int,
)
	R = base_ring(G)
	G.isGB = true
	w = currweight
	ordAlt = G.ord
	while true
		t = next_weight_fractal_walk(G, w, pTargetWeights[p])

		# Handling the final step in the current depth.
		# next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
		# -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.        if t == 1 && p != 1
		if t == 1 && p != 1
			if same_cone(G, T)
				if infoLevel >= 1
					println("depth $p: in cone ", currweight, ".")
				end

				# Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
				if !inCone(G, T, pTargetWeights, p)
					global pTargetWeights =
						[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
					if infoLevel >= 1
						println("depth $p: not in cone ", pTargetWeights[p], ".")
					end
				end
				return G
			end
		elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
			if inCone(G, T, pTargetWeights, p)
				if infoLevel >= 1
					println("depth $p: in cone ", pTargetWeights[p], ".")
				end
				return G
			end
			global pTargetWeights =
				[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
			if infoLevel >= 1
				println("depth $p: not in cone ", pTargetWeights[p], ".")
			end
			continue
		end

		w = w + t * (pTargetWeights[p] - w)
		w = convert_bounding_vector(w)
		Gw = ideal(initials(G, w))

		# Handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
		if !checkInt32(w)
			println("depth $p: Some entries of $w are bigger than Int32. Trying to find another weight,")
			w, b = truncw(G, w, gens(Gw))
			if !b
				println("depth $p: Doing a direct conversion to the target monomial ordering.")
				ordNew = matrix_ordering(R, T)
				w = T[1, :]
				G = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

				if !inCone(G, T, pTargetWeights, p)
					global pTargetWeights =
						[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
					if infoLevel >= 1
						println("depth $p: not in cone ", pTargetWeights[p], ".")
					end
				end
				return G
			end
		end
		ordNew = create_order(R, w, T)
		# Converting the Groebner basis
		if p == nvars(R)
			H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

			if infoLevel >= 1
				println("depth $p: conversion in ", w, ".")
			end
			raise_step_counter()
		else
			if infoLevel >= 1
				println("depth $p: recursive call in $w.")
			end
			H = fractal_recursiv(
				Oscar.IdealGens(R, gens(Gw), ordAlt),
				S,
				T,
				deepcopy(currweight),
				pTargetWeights,
				p + 1,
			)
		end
		#H = liftGW2(G, R, Gw, H, Rn)
		H = lift_fractal_walk(G, H, ordNew)
		G = interreduce_walk(H)
		ordAlt = ordNew
		currweight = w
	end
end

###############################################################
# Extends the plain Fractal Walk by checking the start order.
###############################################################

function fractal_walk_start_order(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
)
	if infoLevel >= 1
		println("fractal_walk_withStartorder results")
		println("Crossed Cones in: ")
	end

	global pTargetWeights =
		[perturbed_vector(G, T, i) for i ∈ 1:nvars(Singular.base_ring(G))]
	return fractal_walk_recursiv_startorder(
		G,
		S,
		T,
		S[1, :],
		pTargetWeights,
		1,
	)
end

function fractal_walk_recursiv_startorder(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
	currweight::Vector{Int},
	pTargetWeights::Vector{Vector{Int}},
	p::Int,
)
	R = base_ring(G)
	G.isGB = true
	ordAlt = G.ord

	# Handling the starting weight.
	if (p == 1)
		if !ismonomial(initials(G, currweight))
			global pStartWeights = [perturbed_vector(G, S, i) for i ∈ 1:nvars(R)]
			global firstStepMode = true
		end
	end
	if firstStepMode
		w = pStartWeights[p]
	else
		w = currweight
	end

	while true
		t = next_weight_fractal_walk(G, w, pTargetWeights[p])

		# Handling the final step in the current depth.
		# next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
		# -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.        if t == 1 && p != 1
		if t == 1 && p != 1
			if same_cone(G, T)
				if infoLevel >= 1
					println("depth $p: in cone ", currweight, ".")
				end

				# Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
				if !inCone(G, T, pTargetWeights, p)
					global pTargetWeights =
						[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
					if infoLevel >= 1
						println("depth $p: not in cone ", pTargetWeights[p], ".")
					end
				end
				return G
			end
		elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
			if inCone(G, T, pTargetWeights, p)
				if infoLevel >= 1
					println("depth $p: in cone ", pTargetWeights[p], ".")
				end
				return G
			end
			global pTargetWeights =
				[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
			if infoLevel >= 1
				println("depth $p: not in cone ", pTargetWeights[p], ".")
			end
			continue
		end

		w = w + t * (pTargetWeights[p] - w)
		w = convert_bounding_vector(w)
		Gw = ideal(initials(G, w))

		# Handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
		if !checkInt32(w)
			w, b = truncw(G, w, gens(Gw))
			if !b
				ordNew = matrix_ordering(R, T)
				w = T[1, :]
				G = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

				if !inCone(G, T, pTargetWeights, p)
					global pTargetWeights =
						[perturbed_vector(G, T, i) for i ∈ 1:nvars(R)]
					if infoLevel >= 1
						println("depth $p: not in cone ", pTargetWeights[p], ".")
					end
				end
				return G
			end
		end
		ordNew = create_order(R, w, T)

		# Converting the Groebner basis
		if p == nvars(R)
			H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

			if infoLevel >= 1
				println("depth $p: conversion in ", w, ".")
			end
			raise_step_counter()
		else
			if infoLevel >= 1
				println("depth $p: recursive call in $w.")
			end
			H = fractal_walk_recursiv_startorder(
				Oscar.IdealGens(R, gens(Gw), ordAlt),
				S,
				T,
				deepcopy(currweight),
				pTargetWeights,
				p + 1,
			)
			global firstStepMode = false
		end
		#H = liftGW2(G, R, Gw, H, Rn)
		H = lift_fractal_walk(G, H, ordNew)
		G = interreduce_walk(H)
		ordAlt = ordNew
		currweight = w
	end
end


###############################################################
# Tran´s version of the Groebner Walk.
# Returns the intermediate Groebner basis if an entry of an intermediate weight vector is bigger than int32.
###############################################################

function tran_walk(
	G::Oscar.IdealGens,
	S::Matrix{Int},
	T::Matrix{Int},
)
	if infoLevel >= 1
		println("tran_walk results")
		println("Crossed Cones in: ")
	end
	currweight = S[1, :]
	tarweight = T[1, :]
	R = base_ring(G)
	if !ismonomial(initials(G, currweight))
		currweight = perturbed_vector(G, S, nvars(R))
	end

	while true
		w = next_weight(G, currweight, tarweight)

		# return the Groebner basis if an entry of w is bigger than int32.
		if !checkInt32(w)
			w, b = truncw(G, w, initials(G, w))
			!b && throw(error("Some entries of the intermediate weight-vector $w are bigger than int32. Choose an other algorithm instead."))
		end
		if w == tarweight
			if same_cone(G, T)
				if infoLevel >= 1
					println("Cones crossed: ", counter)
				end
				return G
			elseif inSeveralCones(initials(G, tarweight))
				tarweight = representation_vector(G, T)
				continue
			end
		end
		G = standard_step_without_int32_check(G, w, T)
		if infoLevel >= 1
			println(w)
			if infoLevel == 2
				println(G)
			end
		end
		currweight = w
		raise_step_counter()
	end
end

###############################################################
# Standard step without checking of the entries of a given weight vector.
###############################################################

function standard_step_without_int32_check(
	G::Oscar.IdealGens,
	w::Vector{Int},
	T::Matrix{Int},
)
	R = base_ring(G)
	ordAlt = G.ord
	ordNew = create_order(R, w, T)
	Gw = initials(G, w)
	H = groebner_basis(ideal(Gw), ordering = ordNew, complete_reduction = true, algorithm = :buchberger)
	#H = liftGW2(G, R, Gw, H, Rn)
	H = lift(G, ordAlt, H, ordNew)
	return interreduce_walk(H)
end


#################################################################
# Procedures of the fractal walk.
# The fractal walk is proposed by Amrhein & Gloor (1998).
#################################################################

# lifts the Groebner basis G to the Groebner basis w.r.t. in the Fractal Walk like it´s done in Fukuda et. al (2005).
function lift_fractal_walk(
	G::Oscar.IdealGens,
	H::Oscar.IdealGens,
	orderingNew::MonomialOrdering,
)
	R = base_ring(G)
	G.isGB = true
	G = Oscar.IdealGens(
		R,
		[
			gen -
			Oscar.IdealGens([reduce(gen, gens(G), ordering = G.ord, complete_reduction = true)], H.ord)[1] for
			gen in gens(H)
		],
		orderingNew,
	)
	G.isGB = true
	return G
end

# returns ´true´ if all polynomials of the given array are monomials.
function ismonomial(Gw::Vector{T}) where {T <: MPolyRingElem}
	for g in Gw
		if length(coefficients(g)) > 1
			return false
		end
	end
	return true
end

# returns ´true´ if all polynomials of the given array are binomials or less.
function isbinomial(Gw::Vector{T}) where {T <: MPolyRingElem}
	for g in Gw
		if length(coefficients(g)) > 2
			return false
		end
	end
	return true
end

# returns the next t to compute the next weight vector w(t) = w + t * (tw - w) like it´s done in Amrhein & Gloor (1998). This Method is NOT tested sufficiently.
function nextT(
	G::Oscar.IdealGens,
	w::Array{T, 1},
	tw::Array{K, 1},
) where {T <: Number, K <: Number}
	if (w == tw)
		return [0]
	end
	tmin = 2
	t = 0
	for g in gens(G)
		a = Singular.leading_exponent_vector(g)
		d = Singular.exponent_vectors(tail(g))
		for v in d
			frac = (dot(w, a) - dot(w, v) + dot(tw, v) - dot(tw, a))
			if frac != 0
				t = (dot(w, a) - dot(w, v)) // frac
			end
			if t > 0 && t < tmin
				tmin = t
			end
		end
	end
	if tmin <= 1
		return tmin
	else
		return [0]
	end
end

# returns the next t to compute the next weight vector w(t) = w + t * (tw - w) like it´s done in Cox, Little & O'Sheao (2005).
function next_weight_fractal_walk(
	G::Oscar.IdealGens,
	cweight::Array{T, 1},
	tweight::Array{K, 1},
) where {T <: Number, K <: Number}
	if (cweight == tweight)
		return [0]
	end
	tmin = 1
	for v in difference_lead_tail(G)
		tw = dot(tweight, v)
		if tw < 0
			cw = dot(cweight, v)
			t = cw // (cw - tw)
			if t < tmin
				tmin = t
			end
		end
	end

	# BigInt is needed to prevent overflows in the conversion of the weight vectors.
	return BigInt(numerator(tmin)) // BigInt(denominator(tmin))
end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G w.r.t the weighted monomial ordering with weight vector of pvecs[p] (pvecs[p-1]) and the matrixordering T.
function inCone(
	G::Oscar.IdealGens,
	T::Matrix{Int},
	pvecs::Vector{Vector{Int}},
	p::Int,
)
	# pvecs(1) equals the first row of T
	if p == 1
		return true
	end
	ord = matrix_ordering(base_ring(G), T)
	cvzip = zip(
		gens(G),
		initials(G, pvecs[p-1]),
		initials(G, pvecs[p]),
	)
	for (g, in, in2) in cvzip
		if !isequal(
			leading_term(g, ordering = ord),
			leading_term(in, ordering = ord),
		) ||
		   !isequal(
			leading_term(g, ordering = ord),
			leading_term(in2, ordering = ord),
		)
			return false
		end
	end
	return true
end


#################################################################
# Procedures of the generic walk.
# The generic walk is proposed by Fukuda, Lauritzen & Thomas (2005).
#################################################################

# returns the initials of the polynomials w.r.t. the vector v.
function facet_initials(
	G::Oscar.IdealGens,
	lm::Vector{T},
	v::Vector{Int},
) where {T <: MPolyRingElem}
	Rn = parent(first(G))
	initials = Array{Singular.elem_type(Rn), 1}(undef, 0)
	count = 1
	for g in G
		inw = Singular.MPolyBuildCtx(Rn)
		el = first(Singular.exponent_vectors(lm[count]))
		for (e, c) in
			zip(Singular.exponent_vectors(g), Singular.coefficients(g))
			if el == e || isparallel(el - e, v)
				Singular.push_term!(inw, c, e)
			end
		end
		h = finish(inw)
		push!(initials, h)
		count += 1
	end
	return initials
end

# returns the differences of the exponent vectors of the leading terms and the polynomials of the generators of I.
function difference_lead_tail(
	G::Oscar.IdealGens,
	Lm::Vector{L},
	T::Union{Matrix{N}, MatElem{N}},
) where {L <: MPolyRingElem, N}
	v = Vector{Int}[]
	for i ∈ 1:length(G)
		ltu = collect(AbstractAlgebra.exponent_vectors(leading_term(Lm[i], ordering = matrix_ordering(base_ring(G), T))))[1]
		for e in AbstractAlgebra.exponent_vectors(G[i])
			if ltu != e
				push!(v, ltu .- e)
			end
		end
	end
	return unique!(v)
end

# returns true if the vector u is parallel to the vector v.
function isparallel(u::Vector{Int}, v::Vector{Int})
	count = 1
	x = 0
	for i ∈ 1:length(u)
		if u[i] == 0
			if v[count] == 0
				count += +1
			else
				return false
			end
		else
			x = v[count] // u[i]
			count += 1
			break
		end
	end
	if count > length(v)
		return true
	end
	for i ∈ count:length(v)
		@inbounds if v[i] != x * u[i]
			return false
		end
	end
	return true
end

# performs the lifting in the generic Walk like it´s proposed by Fukuda et al. (2005).
function lift_generic(
	G::Vector{T},
	Lm::Vector{T},
	H::Vector{T},
	ord::MonomialOrdering,
) where {T <: MPolyRingElem}
	R = parent(first(G))
	Newlm = Array{elem_type(R), 1}(undef, 0)
	liftPolys = Array{elem_type(R), 1}(undef, 0)
	for g in H
		push!(Newlm, leading_term(g, ordering = ord))
		push!(liftPolys, g - reduce_walk(g, G, Lm, ord))
	end
	return liftPolys, Newlm
end

# returns all v \in V if v<0 w.r.t. the ordering represented by T and v>0 w.r.t the ordering represented by S.
function filter_by_ordering(S::Matrix{Int}, T::Matrix{Int}, V::Vector{Vector{Int}})
	btz = Set{Vector{Int}}()
	for v in V
		if less_than_zero(T, v) && bigger_than_zero(S, v)
			push!(btz, v)
		end
	end
	return btz
end

# returns all v \in V if w<v w.r.t. the facet-preorder.
function filter_lf(
	w::Vector{Int},
	S::Matrix{Int},
	T::Matrix{Int},
	V::Set{Vector{Int}},
)
	btz = Set{Vector{Int}}()
	for v in V
		if less_facet(w, v, S, T)
			push!(btz, v)
		end
	end
	return btz
end

# computes the next vector in the generic walk.
function next_gamma(
	G::Oscar.IdealGens,
	Lm::Vector{L},
	w::Vector{Int},
	S::Matrix{Int},
	T::Matrix{Int},
) where {L <: MPolyRingElem}
	V = filter_by_ordering(S, T, difference_lead_tail(G, Lm, T))
	if (w != [0])
		V = filter_lf(w, S, T, V)
	end
	if isempty(V)
		return V
	end
	minV = first(V)
	for v in V
		if less_facet(v, minV, S, T)
			minV = v
		end
	end
	return minV
end

# tests if v>0 w.r.t. the ordering M.
function bigger_than_zero(M::Matrix{Int}, v::Vector{Int})
	nrows, ncols = size(M)
	for i ∈ 1:nrows
		d = 0
		for j ∈ 1:ncols
			@inbounds d += M[i, j] * v[j]
		end
		if d != 0
			return d > 0
		end
	end
	return false
end

# tests if v<0 w.r.t. the ordering M.
function less_than_zero(M::Matrix{Int}, v::Vector{Int})
	nrows, ncols = size(M)
	for i ∈ 1:nrows
		d = 0
		for j ∈ 1:ncols
			@inbounds d += M[i, j] * v[j]
		end
		if d != 0
			return d < 0
		end
	end
	return false
end

# tests if u<v w.r.t. the facet-preorder represented by the matrices S and T.
function less_facet(
	u::Vector{Int},
	v::Vector{Int},
	S::Matrix{Int},
	T::Matrix{Int},
)
	for i ∈ 1:size(T, 1)
		for j ∈ 1:size(S, 1)
			@inbounds Tuv = dot(T[i, :], u) * dot(S[j, :], v)
			@inbounds Tvu = dot(T[i, :], v) * dot(S[j, :], u)
			if Tuv != Tvu
				return Tuv < Tvu
			end
		end
	end
	return false
end

# returns the multiple m for all terms q in p with lm * m = q.
function divides_walk(
	p::MPolyRingElem,
	lm::MPolyRingElem,
	S::MPolyRing,
)
	div = false
	newpoly = MPolyBuildCtx(S)
	for term in terms(p)
		(b, c) = divides(term, lm)
		if b
			push_term!(
				newpoly,
				first(coefficients(c)),
				first(AbstractAlgebra.exponent_vectors(c)),
			)
			div = true
		end
	end
	return finish(newpoly), div
end

# returns p reduced by the Groebner basis G w.r.t. the leading terms Lm.
function reduce_walk(
	p::MPolyRingElem,
	G::Vector{T},
	Lm::Vector{T},
	ord::MonomialOrdering,
) where {T <: MPolyRingElem}
	for i ∈ 1:length(G)
		(q, b) = divides_walk(p, Lm[i], parent(p))
		if b
			return reduce_walk(p - (q * G[i]), G, Lm, ord)
		end
	end
	return p
end

# this function interreduces the Groebner basis G w.r.t. the leading terms Lm with tail-reduction.
function interreduce(
	G::Vector{T},
	Lm::Vector{T},
	ord::MonomialOrdering,
) where {T <: MPolyRingElem}
	for i ∈ 1:Singular.length(G)
		G[i] = reduce_walk(G[i], G[1:end.!=i], Lm[1:end.!=i], ord)
	end
	return G
end

###############################################################
# Several Procedures for the Groebner Walk
###############################################################

# returns the next intermediate weight vector.
function next_weight(
	G::Oscar.IdealGens,
	cw::Vector{Int},
	tw::Vector{Int},
)
	tmin = 1
	for v in difference_lead_tail(G)
		tdotw = dot(tw, v)
		if tdotw < 0
			cdotw = dot(cw, v)
			t = cdotw // (cdotw - tdotw)
			if t < tmin
				tmin = t
			end
		end
	end
	return convert_bounding_vector(cw + BigInt(numerator(tmin)) // BigInt(denominator(tmin)) * (tw - cw))
end

# multiplies every entry of the given weight w with 0.1 as long as it stays on the same halfspace as w.
function truncw(
	G::Oscar.IdealGens,
	w::Vector{Int},
	inw::Vector{T},
) where {T <: MPolyRingElem}
	while !checkInt32(w)
		for i ∈ 1:length(w)
			w[i] = round(w[i] * 0.10)
		end
		w = convert_bounding_vector(w)
		if inw != initials(G, w)
			# initials are different - return unrounded weight
			return w, false
		end
	end
	# converted to Vector w of the same face
	return w, true
end

# returns the initialform of G w.r.t. the given weight vector.
function initials(
	G::Oscar.IdealGens,
	w::Vector{Int},
)
	R = base_ring(G)
	inits = elem_type(R)[]
	for g in gens(G)
		inw = MPolyBuildCtx(R)
		maxw = 0
		eczip = zip(AbstractAlgebra.exponent_vectors(g), coefficients(g))
		for (e, c) in eczip
			tmpw = dot(w, e)
			if maxw == tmpw
				push_term!(inw, c, e)
			elseif maxw < tmpw
				inw = MPolyBuildCtx(R)
				push_term!(inw, c, e)
				maxw = tmpw
			end
		end
		h = finish(inw)
		push!(inits, h)
	end
	return inits
end

function difference_lead_tail(I::Oscar.IdealGens)
	v = Vector{Int}[]
	for g in gens(I)
		leadexpv = first(AbstractAlgebra.exponent_vectors(leading_term(g, ordering = I.ord)))
		tailexpvs = AbstractAlgebra.exponent_vectors(tail(g, ordering = I.ord))
		for e in tailexpvs
			push!(v, leadexpv .- e)
		end
	end
	return unique!(v)
end

# computes a p-perturbed vector from the matrix M.
function perturbed_vector(G::Oscar.IdealGens, M::Matrix{Int}, p::Integer)
	m = Int[]
	n = size(M, 1)
	for i ∈ 1:p
		max = M[i, 1]
		for j ∈ 1:n
			temp = abs(M[i, j])
			if temp > max
				max = temp
			end
		end
		push!(m, max)
	end
	msum = 0
	for i ∈ 2:p
		msum += m[i]
	end
	maxdeg = 0
	for g in gens(G)
		td = deg(g, n)
		if (td > maxdeg)
			maxdeg = td
		end
	end
	e = maxdeg * msum + 1
	w = M[1, :] * e^(p - 1)
	for i ∈ 2:p
		w += e^(p - i) * M[i, :]
	end
	return convert_bounding_vector(w)
end

# returns 'true' if the leading terms of G w.r.t the matrixorder T are the same as the leading terms of G w.r.t the weighted monomial order with weight vector t and matrix T.
#function inCone(G::Oscar.IdealGens, T::Union{Matrix{N}, MatElem{N}}, t::Vector{Int})
#    R = change_order(G.base_ring, T)
#    I = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
#    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), t))
#    for (g, ing) in cvzip
#        if !isequal(Singular.leading_exponent_vector(g), Singular.leading_exponent_vector(ing))
#            return false
#        end
#    end
#    return true
#end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G w.r.t the weighted monomial order with weight vector t and the matrix order T.
#function inCone(G::Oscar.IdealGens, t::Vector{Int})
#    cvzip = zip(Singular.gens(G), initials(base_ring(G), gens(G), t))
#    for (g, ing) in cvzip
#        if !isequal(leading_exponent_vector(g), leading_exponent_vector(ing))
#            return false
#        end
#    end
#    return true
#end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G with the current ordering.
function same_cone(G::Oscar.IdealGens, T::Matrix{Int})
	R = base_ring(G)
	ord = matrix_ordering(R, T)
	for g in gens(G)
		if leading_term(g, ordering = ord) != leading_term(g, ordering = G.ord)
			return false
		end
	end
	return true
end

# lifts the Groebner basis G to the Groebner basis w.r.t. the Ring Rn like it´s presented in Fukuda et al. (2005).
function lift(
	G::Oscar.IdealGens,
	orderingAlt::MonomialOrdering,
	H::Oscar.IdealGens,
	ordering::MonomialOrdering,
)
	G = Oscar.IdealGens(
		[
			gen - Oscar.IdealGens([reduce(gen, gens(G), ordering = orderingAlt, complete_reduction = true)], ordering)[1]
			for gen in gens(H)
		],
		ordering,
		isGB = true,
	)
	return G
end

# lifts the Groebner basis G to the Groebner basis w.r.t. the Ring Rn like it´s done in Collart et al. (1997).
# FIXME: Needs to be fixed.
function liftGW2(
	G::Oscar.IdealGens,
	orderingAlt::MonomialOrdering,
	inG::MPolyIdeal{T},
	H::Oscar.IdealGens,
	ordering::MonomialOrdering,
) where {T <: MPolyRingElem}
	H = gens(H)
	R = base_ring(G)
	G = gens(G)
	for i ∈ 1:length(H)
		q = divrem(gens(Oscar.IdealGens([H[i]], orderingAlt)), gens(inG))
		H[i] = R(0)
		for j ∈ 1:length(gens(inG))
			println(gens(inG)[j])
			println(G[j])
			H[i] = H[i] + q[1][1][j] * G[j]
		end
	end
	return Oscar.IdealGens(filter(!iszero, H), ordering)
end

# divisionalgorithm that returns q with q_1*f_1 + ... + q_s *f_s=p.
function division_algorithm(
	p::T,
	f::Vector{T},
	R::MPolyRing,
) where {T <: MPolyRingElem}
	q = Array{Singular.elem_type(R), 1}(undef, length(f))
	for i ∈ 1:length(f)
		q[i] = R(0)
	end
	while !isequal(p, R(0))
		i = 1
		div = false
		while (div == false && i <= length(f))
			b, m = divides(leading_term(p), leading_term(f[i]))
			if b
				q[i] = q[i] + m
				p = p - (m * f[i])
				div = true
			else
				i = i + 1
			end
		end
		if div == false
			p = p - leading_term(p)
		end
	end
	return q
end

# converts a vector wtemp by dividing the entries with gcd(wtemp).
function convert_bounding_vector(wtemp::Vector{T}) where {T <: Rational{BigInt}}
	w = Vector{Int}()
	g = gcd(wtemp)
	for i ∈ 1:length(wtemp)
		push!(w, float(divexact(wtemp[i], g)))
	end
	return w
end

# converts a vector wtemp by dividing the entries with gcd(wtemp).
function convert_bounding_vector(wtemp::Vector{T}) where {T <: Number}
	w = Vector{Int}()
	g = gcd(wtemp)
	for i ∈ 1:length(wtemp)
		push!(w, float(divexact(wtemp[i], g)))
	end
	return w
end

# returns a copy of the PolynomialRing I, equipped with the ordering weight_ordering(cw)*matrix_ordering(T).
function create_order(
	R::MPolyRing,
	cw::Array{L, 1},
	T::Matrix{Int},
) where {L <: Number}
	s = size(T)
	if s[1] == s[2]
		ord = weight_ordering(cw, matrix_ordering(R, T))
	elseif s[1] - s[2] == 1
		ord = weight_ordering(cw,
			matrix_ordering(R, T))
	else
		ord = weight_ordering(cw,
			matrix_ordering(R, T),
		)
	end
	return ord
end

# interreduces the Groebner basis G.
function interreduce_walk(G::Oscar.IdealGens)
	Rn = base_ring(G)
	Generator = collect(gens(G))
	I = 0
	for i ∈ 1:length(gens(G))
		Generator[i] = reduce(Generator[i], Generator[1:end.!=i], ordering = G.ord, complete_reduction = true)
	end
	return Oscar.IdealGens(Generator, G.ord)
end

#############################################
# unspecific help functions
#############################################

function ident_matrix(n::Int)
	M = zeros(Int, n, n)
	for i ∈ 1:n
		M[i, i] = 1
	end
	return M
end

function anti_diagonal_matrix(n::Int)
	M = zeros(Int, n, n)
	for i ∈ 1:n
		M[i, n+1-i] = -1
	end
	return M
end

# Singular.isequal depends on order of generators
function equalitytest(G::Oscar.IdealGens, K::Oscar.IdealGens)
	if length(gens(G)) != length(gens(K))
		return false
	end
	generators = Singular.gens(G)
	count = 0
	for gen in generators
		for r in Singular.gens(K)
			if gen * first(coefficients(leading_term(r, ordering = G.ord))) - r * first(coefficients(leading_term(gen, ordering = G.ord))) == 0
				count += 1
				break
			end
		end
	end
	if count == length(gens(G))
		return true
	end
	return false
end

function dot(v::Vector{Int}, w::Vector{Int})
	sum = 0
	for i ∈ 1:length(v)
		@inbounds sum += v[i] * w[i]
	end
	return sum
end

function ordering_as_matrix(w::Vector{Int}, ord::Symbol)
	if length(w) > 2
		if ord == :lex
			return [
				w'
				ident_matrix(length(w))[1:length(w)-1, :]
			]
		end
		if ord == :deglex
			return [
				w'
				ones(Int, length(w))'
				ident_matrix(length(w))[1:length(w)-2, :]
			]
		end
		if ord == :degrevlex
			return [
				w'
				ones(Int, length(w))'
				anti_diagonal_matrix(length(w))[1:length(w)-2, :]
			]
		end
		if ord == :revlex
			return [
				w'
				anti_diagonal_matrix(length(w))[1:length(w)-1, :]
			]
		end
	else
		error("not implemented")
	end
end

function change_weight_vector(w::Vector{Int}, M::Matrix{Int})
	return [
		w'
		M[2:length(w), :]
	]
end

function insert_weight_vector(w::Vector{Int}, M::Matrix{Int})
	return [
		w'
		M[1:length(w)-1, :]
	]
end

function add_weight_vector(w::Vector{Int}, M::Matrix{Int})
	return [
		w'
		M
	]
end

function ordering_as_matrix(ord::Symbol, nvars::Int)
	if ord == :lex
		return ident_matrix(nvars)
	end
	if ord == :deglex
		return [
			ones(Int, nvars)'
			ident_matrix(nvars)[1:nvars-1, :]
		]
	end
	if ord == :degrevlex
		return [
			ones(Int, nvars)'
			anti_diagonal_matrix(nvars)[1:nvars-1, :]
		]
	end
	if ord == :revlex
		return [
			w'
			anti_diagonal_matrix(length(w))[1:length(w)-1, :]
		]
	else
		error("not implemented")
	end
end

function deg(p::MPolyRingElem, n::Int)
	max = 0
	for mon in Singular.monomials(p)
		ev = Singular.exponent_vectors(mon)
		sum = 0
		for e in ev
			for i ∈ 1:n
				sum += e[i]
			end
			if (max < sum)
				max = sum
			end
		end
	end
	return max
end


function checkInt32(w::Vector{Int})
	for i ∈ 1:length(w)
		if tryparse(Int32, string(w[i])) == nothing
			return false
		end
	end
	return true
end

# computes the representation of the matrixorder defined by T.
function representation_vector(G::Oscar.IdealGens, T::Matrix{Int})
	n = size(T)[1]
	M = 0
	for i ∈ 1:n
		for j ∈ 1:n
			temp = T[i, j]
			if M < temp
				M = temp
			end
		end
	end
	d0 = 0
	for g in Singular.gens(G)
		temp = deg(g, n)
		if d0 < temp
			d0 = temp
		end
	end
	d = M * (2 * d0^2 + (n + 1) * d0)
	w = d^(n - 1) * T[1, :]
	for i ∈ 2:n
		w = w + d^(n - i) * T[i, :]
	end
	return w
end

# checks if Gw is an initialform of an ideal corresponding to a face of the Groebner fan with dimension less than n-1.
function inSeveralCones(Gw::Vector{T}) where {T <: MPolyRingElem}
	counter = 0
	for g in Gw
		if size(collect(Singular.coefficients(g)))[1] > 2
			return true
		end
		if size(collect(Singular.coefficients(g)))[1] == 2
			counter = counter + 1
		end
	end
	if counter > 1
		return true
	end
	return false
end