@testset "Testing the Groebner Walk" begin
	R, (a, b, c, d) = polynomial_ring(
		Singular.N_ZpField(32003),
		["a", "b", "c", "d"],
		ordering = :degrevlex,
	)
	id = ideal(
		R,
		[
			2 * a^2 * b +
			3 * a * b^2 +
			3 * b^3 +
			4 * c^3 +
			4 * a * b * d +
			c^2 * d +
			2 * b * d^2 +
			2 * d^3 +
			4 * c^2 +
			2 * c * d +
			2 * c,
			2 * a^2 * b +
			5 * a^2 * c +
			2 * b * c^2 +
			c^2 * d +
			a * c +
			2 * b * d +
			2 * c * d +
			d^2 +
			a +
			4 * d,
		],
	)

	ideals = []
	infoLevel = 1
	for i ∈ 2:nvars(R)
		push!(ideals, groebner_walk(id, degrevlex(R), lex(R), :perturbed, i))
	end
	push!(ideals, groebner_walk(id, degrevlex(R), lex(R), :standard))
	# push!(ideals, groebnerwalk(id, degrevlex(R), lex(R), :tran))

	push!(ideals, groebnerwalk(id, degrevlex(R), lex(R), :fractal))
	push!(ideals, groebnerwalk(id, degrevlex(R), lex(R), :fractal_start_order))
	push!(ideals, groebnerwalk(id, degrevlex(R), lex(R), :fractal_combined))
	push!(ideals, groebnerwalk(id, degrevlex(R), lex(R), :generic))

	s = groebner_basis(id, ordering = lex(R), complete_reduction = true)
	for i in ideals
		i.isGB = true
		@test equalitytest(
			Oscar.IdealGens(R, gens(i), lex(R)),
			s,
		)
	end

	R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
	I = ideal([y^4 + x^3 - x^2 + x, x^4])
	id = groebnerwalk(I, degrevlex(R), lex(R), :tran)

	s = groebner_basis(I, ordering = lex(R), complete_reduction = true)
	@test equalitytest(
		Oscar.IdealGens(R, gens(id), lex(R)),
		s,
	)

	@testset "backend-functions for the Groebner Walk" begin

		R, (x, y, z) = polynomial_ring(
			QQ,
			["x", "y", "z"],
			ordering = :degrevlex,
		)

		f1 = 3 * x^2 + 16 * x^2 * z + 14 * x^2 * y^3
		f2 = y^3 * z + 17 * x^2 * z^2 + 7 * x^2 * y^2 * z^2 + 13 * x^3 * z^2
		I = Oscar.IdealGens(R, [f1, f2], degrevlex(R))
		sol = [14 * x^2 * y^3, y^3 * z + 7 * x^2 * y^2 * z^2]
		@test initials(I, [1, 3, 1]) == sol

		@test difference_lead_tail(I) ==
			  [[0, 3, -1], [0, 3, 0], [-1, 2, 0], [2, -1, 1], [0, 2, 0]]

		F = [
			13 * x^3 * z^2,
			14 * x^2 * y^3,
			98 * x * y^5 * z^2,
			y^7 * z + x^2 * z^3,
			14 * x * y^6 * z,
		]
		g = y^7 * z + x^2 * z^3 + 28 * x^2 * y^4
		q = Array{elem_type(R), 1}(undef, 5)
		q[1] = R(0)
		q[2] = R(2 * y)
		q[3] = R(0)
		q[4] = R(1)
		q[5] = R(0)
		@test division_algorithm(g, F, R) == q

		J = Oscar.IdealGens(R, [f2, f1])
		f1 = 4 * x^2 + 16 * x^2 * z + 14 * x^2 * y^3
		f2 = y^3 * z + 17 * x^2 * z^2 + 7 * x^2 * y^2 * z^2 + 13 * x^3 * z^2
		K = Oscar.IdealGens(R, [f1, f2])
		@test equalitytest(I, J) == true
		@test equalitytest(I, K) == false

		@test deg(f1, 3) == 5

		R, (x, y, z, t, u, v) = polynomial_ring(
			Singular.N_ZpField(32003),
			["x", "y", "z", "t", "u", "v"],
			ordering = :degrevlex,
		)
		StartOrd = ordering_as_matrix(:degrevlex, 6)
		TarOrd = ordering_as_matrix(:lex, 6)
		f1 = 45 * y + 35 * u - 165 * v - 36
		f2 = 36 * y + 25 * z + 40 * t - 27 * u
		f3 = 25 * y * u - 165 * v^2 + 15 * x - 18 * z + 30t
		f4 = 15 * y * z + 20 * t * u - 9 * x
		f5 = -11 * v^3 + x * y + 2 * z * t
		f6 = -11 * u * v + 3 * v^2 + 99 * x
		id = ideal(R, [f1, f2, f3, f4, f5, f6])
		I = groebner_basis(id, complete_reduction = true)
		@test perturbed_vector(I, StartOrd, 1) == [1, 1, 1, 1, 1, 1]
		@test perturbed_vector(I, StartOrd, 2) == [4, 4, 4, 4, 4, 3]
		@test perturbed_vector(I, StartOrd, 3) == [49, 49, 49, 49, 48, 42]
		@test perturbed_vector(I, StartOrd, 4) ==
			  [1000, 1000, 1000, 999, 990, 900]
		@test perturbed_vector(I, TarOrd, 1) == [1, 0, 0, 0, 0, 0]
		@test perturbed_vector(I, TarOrd, 2) == [4, 1, 0, 0, 0, 0]
		@test perturbed_vector(I, TarOrd, 3) == [49, 7, 1, 0, 0, 0]
		@test perturbed_vector(I, TarOrd, 4) == [1000, 100, 10, 1, 0, 0]

		@test same_cone(I, add_weight_vector([1000, 1000, 1000, 999, 990, 900], TarOrd)) == true
		@test same_cone(I, add_weight_vector([100, 1000, 1000, 999, 990, 900], TarOrd)) == false

		@test dot([1, 2, 3, 4], [2, 2, 2, 2]) == 20

		R, (x, y, z) = polynomial_ring(
			QQ,
			["x", "y", "z"],
			ordering = :degrevlex,
		)

		f1 = 3 * x^2
		f2 = y^3 * z + 17 * x^2 * z^2
		I = ideal(R, [f1, f2])
		f1 = 3 * x^2
		f2 = y^3 * z
		J = ideal(R, [f1, f2])
		f1 = x^3 * y^2 + z^2 + y^2
		f2 = y^3
		K = ideal(R, [f1, f2])

		@test ismonomial(gens(I)) == false
		@test ismonomial(gens(J)) == true
		@test isbinomial(gens(I)) == true
		@test isbinomial(gens(J)) == true
		@test isbinomial(gens(K)) == false
		@test ismonomial(gens(K)) == false

		R, (x, y) = polynomial_ring(
			QQ,
			["x", "y"],
			ordering = :degrevlex,
		)

		f1 = x^2 - y^3
		f2 = x^3 - y^2 - x
		I = ideal(R, [f1, f2])
		G = groebner_basis(I, complete_reduction = true)

		@test next_gamma(
			G,
			[leading_term(g) for g in gens(G)],
			[0],
			ordering_as_matrix(:degrevlex, 2),
			ordering_as_matrix(:lex, 2),
		) == [-2, 3]


		g = [R(y^3 - x^2), R(x^3)]
		G.ord = lex(R)
		@test facet_initials(
			G,
			[leading_term(g) for g in G],
			[-2, 3],
		) == g
		G.ord = degrevlex(R)
		@test difference_lead_tail(
			G,
			[leading_term(g) for g in gens(G)], ordering_as_matrix(:lex, 2),
		) == [[-2, 3], [3, -2], [2, 0]]

		@test isparallel([1, 2], [1, 4]) == false
		@test isparallel([1, 2], [2, 4]) == true
		@test isparallel([-1, 0], [-2, 1]) == false
		@test isparallel([-1, 0], [2, 0]) == true

		@test less_facet(
			[-2, 3],
			[-1, 4],
			ordering_as_matrix(:degrevlex, 2),
			ordering_as_matrix(:lex, 2),
		) == true
		@test less_facet(
			[-1, 7],
			[-1, 4],
			ordering_as_matrix(:degrevlex, 2),
			ordering_as_matrix(:lex, 2),
		) == false

		R, (a, b, c, d, e) = PolynomialRing(
			QQ,
			["a", "b", "c", "d", "e"],
			ordering = :degrevlex,
		)
		J = ideal(
			R,
			[
				b + 3 * b^3 + 2 * b * c * e + 5 * b * d * e,
				4 + b^2 + 4 * b * c + 5 * b^3 + c * d * e,
				d * e + 5 * b^2 * e,
			],
		)
		I = groebner_basis(J, complete_reduction = true)
		f1 = R(
			a^3 +
			a^2 +
			b^5 * a^3 * c^9 +
			e^3 +
			b^2 * a^2 * c^4 +
			d^3 +
			e^3 +
			b^2 * d^2 * e^7,
		)
		f2 = R(a^3 * b^3)
		f3 = R(a^2 * b^4 + c^2 + a^3 * 4 + a * e^3)
		f4 = R(0)

		@test (reduce(f3, gens(I), ordering = degrevlex(R), complete_reduction = true)) ==
			  reduce_walk(f3, gens(I), [leading_term(g) for g in gens(I)], degrevlex(R))
		@test (reduce(f2, gens(I), ordering = degrevlex(R), complete_reduction = true)) ==
			  reduce_walk(f2, gens(I), [leading_term(g) for g in gens(I)], degrevlex(R))
		@test (reduce(f2, gens(I), ordering = degrevlex(R), complete_reduction = true)) ==
			  reduce_walk(f2, gens(I), [leading_term(g) for g in gens(I)], degrevlex(R))
		@test (reduce(f4, gens(I), ordering = degrevlex(R), complete_reduction = true)) ==
			  reduce_walk(f4, gens(I), [leading_term(g) for g in gens(I)], degrevlex(R))

		J = groebner_basis(J, ordering = degrevlex(R), complete_reduction = false)
		@test equalitytest(
			groebner_basis(ideal(R, gens(J)), complete_reduction = true), Oscar.IdealGens(interreduce(
				gens(J),
				[leading_term(g) for g in gens(J)],
				degrevlex(R),
			)),
		)
	end
end