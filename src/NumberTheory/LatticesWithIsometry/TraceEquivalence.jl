export trace_lattice, inverse_trace_lattice

@doc Markdown.doc"""
    trace_lattice(L::Hecke.AbsLat{T}; order::Integer = 2) where T -> LatticeWithIsometry

Given a lattice `L` which is either a $\mathbb Z$-lattice or a hermitian lattice
over the $E/K$ with `E` a cyclotomic field and `K` a maximal real subfield of
`E`, return the trace lattice `Lf` whose underlying lattice is the trace lattice
of `L` and the underlying map `f` is:
  - the identity if `L` is a `ZLat` and `order == 1`;
  - the opposite of the identity if `L` is a `ZLat` and `order == 2`;
  - given by the multiplication by a $n$-th root of the unity if `L` is a
    `HermLat`, where `n` is the order a primitive element in `E`.

Note that the optional argument `order` has no effect if `L` is not a
$\mathbb Z$-lattice.
"""
function trace_lattice(L::Hecke.AbsLat{T}; order::Integer = 2) where T
  @req order > 0 "The order must be positive"
  n = degree(L)
  E = base_field(L)
  G = gram_matrix(rational_span(L))

  if E == QQ || E == ZZ
    if order == 1
      f = identity_matrix(E, n)
    elseif order == 2
      f = -identity_matrix(E, n)
    else
      error("For ZLat the order must be 1 or 2")
    end
    return lattice_with_isometry(L, f, order, check = false)
  end
  
  bool, m = Hecke.is_cyclotomic_type(E)

  @req bool "Base field must be cyclotomic"

  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  EtoEabs = pseudo_inv(EabstoE)
  d = degree(Eabs)
  Gabs = map_entries(EtoEabs, G)
  z = gen(Eabs)
  chi = minpoly(z)
  Mchi = companion_matrix(chi)
  f = block_diagonal_matrix([Mchi for i=1:n])

  coeffs = fmpq[]
  for i=1:n
    for iz=1:d
      for j=1:n
        for jz=1:d
	        g = z^(iz-jz)*Gabs[i,j]
	        push!(coeffs, trace(g))
	      end 
      end
    end
  end
  gram = matrix(QQ, d*n, d*n, coeffs)
  @assert f*gram*transpose(f) == gram
  LL = Zlattice(gram = 1//m*gram)
  return lattice_with_isometry(LL, f, m, check = false)
end

@doc Markdown.doc"""
    inverse_trace_lattice(L::ZLat, f::fmpq_mat; n::Integer = -1,
				   check::Bool = true) -> HermLat

Given a $\mathbb Z$-lattice `L` and a matrix with rational entries `f`,
representing an isometry of `L` of order $n \geq 3$, such that the totient of `n`
divides the rank of `L` and the minimal polynomial of `f` is the `n`th cyclotomic
polynomial, return the corresponding hermitian lattice over $E/K$ where `E` is the
`n`th cyclotomic field and `K` is a maximal real subfield of `E`, where the
multiplication by the `n`th root of unity correspond to the mapping via `f`.
"""
function inverse_trace_lattice(L::ZLat, f::fmpq_mat; n::Integer = -1, check::Bool = true,
                                                     ambient_representation::Bool = true)
  if n <= 0
    Oscar._exponent(f)
  end
  
  @req n > 0 "f is not of finite exponent"

  if ambient_representation
    ok, f = can_solve_with_solution(basis_matrix(L), basis_matrix(L)*f, side =:left)
    @req ok "Isometry does not restrict to L"
  end

  if check
    @req n >= 3 "No hermitian inverse trace lattice for order less than 3"
    G = gram_matrix(L)
    @req Oscar._is_cyclotomic_polynomial(minpoly(f)) "The minimal polynomial of f must be irreducible and cyclotomic"
    @req f*G*transpose(f) == G "f does not define an isometry of L"
    @req Oscar._exponent(f) == n "The order of f should be equal to n"
    @req divides(rank(L), euler_phi(n))[1] "The totient of n must divides the rank of L"
  end

  E,b = cyclotomic_field_as_cm_extension(n, cached = false)

  m = divexact(rank(L), euler_phi(n))
  # We choose as basis for the hermitian lattice Lh the identity matrix
  gram = matrix(zeros(E, m,m))
  G = gram_matrix(L)
  v = zero_matrix(QQ, 1, rank(L))

  for i=1:m
    for j=1:m
      vi = deepcopy(v)
      vi[1,1+(i-1)*euler_phi(n)] = one(QQ)
      vj = deepcopy(v)
      vj[1,1+(j-1)*euler_phi(n)] = one(QQ)
      alpha = sum([(vi*G*transpose(vj*f^k))[1]*b^k for k in 0:n-1])
      gram[i,j] = alpha
    end
  end

  s = involution(E)
  @assert transpose(map_entries(s, gram)) == gram

  Lh = hermitian_lattice(E, gram = gram)
  return Lh
end

