export trace_lattice

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
function trace_lattice(L::Hecke.AbsLat{T}; alpha::FieldElem = one(base_field(L)),
                                           beta::FieldElem = gen(base_field(L)),
                                           order::Integer = 2) where T
  E = base_field(L)
  @req maximal_order(E) == equation_order(E) "Equation order and maximal order must coincide"
  @req order > 0 "The order must be positive"
  @req degree(L) == rank(L) "Lattice must be of full rank"
  @req parent(beta) === E "beta must be an element of the base algebra of L"
  n = degree(L)
  s = involution(E)
  if s(beta)*beta != 1
    beta = beta//s(beta)
  end

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
  @req !bool || findfirst(i -> isone(beta^i), 1:m) == m "beta must be a $m-primitive root of 1"

  Lres, f = restrict_scalars_with_map(L, QQ, alpha)
  iso = zero_matrix(QQ, 0, degree(Lres))
  v = vec(zeros(QQ, 1, degree(Lres)))

  for i in 1:degree(Lres)
    v[i] = one(QQ)
    v2 = f(v)
    v2 = beta.*v2
    v3 = f\v2
    iso = vcat(iso, transpose(matrix(v3)))
    v[i] = zero(QQ)
  end

  return lattice_with_isometry(Lres, iso, ambient_representation=true, check = true)
end

function absolute_representation_matrix(b::Hecke.NfRelElem)
  E = parent(b)
  n = absolute_degree(E)
  B = absolute_basis(E)
  m = zero_matrix(QQ, n, n)
  for i in 1:n
    bb = B[i]
    v = absolute_coordinates(b*bb)
    m[i,:] = transpose(matrix(v))
  end
  return m
end

function _hermitian_structure(L::ZLat, f::fmpq_mat; E = nothing,
                                                    n::Integer = -1,
                                                    check::Bool = true,
                                                    ambient_representation::Bool = true)
  if n <= 0
    n = multiplicative_order(f)
  end
  
  @req n > 0 "f is not of finite exponent"

  if ambient_representation
    ok, f = can_solve_with_solution(basis_matrix(L), basis_matrix(L)*f, side =:left)
    @req ok "Isometry does not restrict to L"
  end

  if check
    @req n >= 3 "No hermitian structure for order less than 3"
    G = gram_matrix(L)
    @req is_cyclotomic_polynomial(minpoly(f)) "The minimal polynomial of f must be irreducible and cyclotomic"
    @req f*G*transpose(f) == G "f does not define an isometry of L"
    @req multiplicative_order(f) == n "The order of f should be equal to n"
    @req divides(rank(L), euler_phi(n))[1] "The totient of n must divides the rank of L"
  end
  
  if E === nothing
    E, b = cyclotomic_field_as_cm_extension(n)
  elseif !Hecke.is_cyclotomic_type(E)[1]
    @req degree(E) == 2 && absolute_degree(E) == 2*euler_phi(n) "E should be the $n-th cyclotomic field seen as a cm-extension of its real cyclotomic subfield"
    Et, t = E["t"]
    rt = roots(t^n-1)
    @req length(rt) == euler_phi(n) "E is not of cyclotomic type"
    b = rt[1]
  else
    b = gen(E)
  end
 
  mb = absolute_representation_matrix(b)
  m = divexact(rank(L), euler_phi(n))
  diag = [mb for i in 1:m]
  mb = block_diagonal_matrix(diag)
  bca = Hecke._basis_of_commutator_algebra(f, mb)
  @assert !is_empty(bca)
  l = inv(bca[1])
  B = matrix(absolute_basis(E))
  gene = Vector{elem_type(E)}[]
  for i in 1:nrows(l)
    vv = vec(zeros(E, 1, m))
    v = l[i,:]
    for j in 1:m
      a = (v[1, 1+(j-1)*euler_phi(n):j*euler_phi(n)]*B)[1]
      vv[j] = a
    end
    push!(gene, vv)
  end
  # We choose as basis for the hermitian lattice Lh the identity matrix
  gram = matrix(zeros(E, m, m))
  G = inv(l)*gram_matrix_of_rational_span(L)*inv(transpose(l))
  v = zero_matrix(QQ, 1, rank(L))
  for i=1:m
    for j=1:m
      vi = deepcopy(v)
      vi[1,1+(i-1)*euler_phi(n)] = one(QQ)
      vj = deepcopy(v)
      vj[1,1+(j-1)*euler_phi(n)] = one(QQ)
      alpha = sum([(vi*G*transpose(vj*mb^k))[1]*b^k for k in 0:n-1])
      gram[i,j] = alpha
    end
  end
  s = involution(E)
  @assert transpose(map_entries(s, gram)) == gram

  Lh = hermitian_lattice(E, gene, gram = 1//n*gram)
  return Lh
end

