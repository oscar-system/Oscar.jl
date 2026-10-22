@doc raw"""
    groebner_basis_hilbert_driven(I::MPolyIdeal{P}; destination_ordering::MonomialOrdering,
                    complete_reduction::Bool = false,
                    weights::Vector{Int} = ones(Int, ngens(base_ring(I))),
                    hilbert_numerator::Union{Nothing, ZZPolyRingElem} = nothing) 
                    where {P <: MPolyRingElem}

Return a Gröbner basis of `I` with respect to `destination_ordering`.

!!! note
    The function implements a version of the Hilbert driven Gröbner basis algorithm.
    See the corresponding section of the OSCAR documentation for some details.

!!! note
    The destination ordering must be global.

!!! note
    All weights must be positive. If no weight vector is entered by the user, all weights 
    are set to 1. An error is thrown if the generators of `I` are not homogeneous with 
    respect to the corresponding (weighted) degree. The Hilbert driven Gröbner basis algorithm
    can also be called via the function `groebner_basis`: use `algorithm = :hilbert`.
    Then, the algorithm first tries to find an appropriate weight vector for the given ideal.
    If no such vector exists, the algorithm proceeds via homogenization and dehomogenization.

!!! note
    If $R$ is the parent ring of $I$, and $p, q\in\mathbb Z[t]$ are polynomials
    such that $p/q$ represents the Hilbert series of $R/I$ as a rational function with 
    denominator $q = (1-t^{w_1})\cdots (1-t^{w_n}),$ where $n$ is the number of variables 
    of $R$, and $w_1, \dots, w_n$ are the assigned weights, then `hilbert_numerator` is 
    meant to be $p$.

!!! warning
    If `hilbert_numerator` is set by the user, the function does NOT check
    whether it is a valid numerator for the given data. If `hilbert_numerator`
    is not set by the user, it will be computed internally.

# Examples
```jldoctest
julia> R, (a, b, c, d, e, f, g) = polynomial_ring(QQ, [:a, :b, :c, :d, :e, :f, :g]);

julia> V = [-3*a^2+2*f*b+3*f*d, (3*g*b+3*g*e)*a-3*f*c*b,
                      -3*g^2*a^2-c*b^2*a-g^2*f*e-g^4, e*a-f*b-d*c];

julia> I = ideal(R, V);

julia> o = degrevlex([a, b, c])*degrevlex([d, e, f, g]);

julia> G = groebner_basis_hilbert_driven(I, destination_ordering = o);

julia> length(G)
296

julia> total_degree(G[49])
30
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(GF(32003), [:x, :y, :z]);

julia> f1 = x^2*y+169*y^21+151*x*y*z^10;

julia> f2 = 6*x^2*y^4+x*z^14+3*z^24;

julia> f3 = 11*x^3+5*x*y^10*z^10+2*y^20*z^10+y^10*z^20;

julia> I = ideal(R, [f1, f2,f3]);

julia> W = [10, 1, 1];

julia> GB = groebner_basis_hilbert_driven(I, destination_ordering = lex(R), weights = W);

julia> length(GB)
40
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(GF(32003), [:x, :y, :z]);

julia> f1 = x^2*y+169*y^21+151*x*y*z^10;

julia> f2 = 6*x^2*y^4+x*z^14+3*z^24;

julia> f3 = 11*x^3+5*x*y^10*z^10+2*y^20*z^10+y^10*z^20;

julia> I = ideal(R, [f1, f2,f3]);

julia> W = [10, 1, 1];

julia> S, t = polynomial_ring(ZZ, :t)
(Univariate polynomial ring in t over ZZ, t)

julia> hn = -t^75 + t^54 + t^51 + t^45 - t^30 - t^24 - t^21 + 1
-t^75 + t^54 + t^51 + t^45 - t^30 - t^24 - t^21 + 1

julia> GB = groebner_basis_hilbert_driven(I, destination_ordering = lex(R), weights = W, hilbert_numerator = hn);

julia> length(GB)
40
```
"""
function groebner_basis_hilbert_driven(I::MPolyIdeal{P};
                                       destination_ordering::MonomialOrdering,
                                       complete_reduction::Bool = false,
                                       weights::Vector{Int} = ones(Int, ngens(base_ring(I))),
                                       hilbert_numerator::Union{Nothing, ZZPolyRingElem} = nothing) where {P <: MPolyRingElem}

  ord_dest = destination_ordering
  haskey(I.gb, ord_dest) && return I.gb[ord_dest]
  
  isa(coefficient_ring(I), AbstractAlgebra.Field) || error("The underlying coefficient ring of I must be a field.")
  is_global(ord_dest) || error("Destination ordering must be global.")
  all(f -> _is_homogeneous(f, weights), gens(I)) || error("I must be given by generators homogeneous with respect to the given weights.")

  R = base_ring(I)
  ord_start = all(isone, weights) ? degrevlex(R) : wdegrevlex(R, weights)
  if isnothing(hilbert_numerator)
    G = standard_basis(I, ordering = ord_start)
    h = Singular.hilbert_series(singular_generators(G, G.ord), weights)
  else
    h = (Int32).([coeff(hilbert_numerator, i) for i in 0:degree(hilbert_numerator)+1])
    #TODO Is this still needed? Should we handle this differently?
  end

  singular_I_gens = singular_generators(I.gens, ord_dest)
  singular_ring = base_ring(singular_I_gens)
  SI = Singular.Ideal(singular_ring, gens(singular_I_gens)...)
  GBS  = Singular.std_hilbert(SI, h, (Int32).(weights),
                            complete_reduction = complete_reduction)
  GB = IdealGens(base_ring(I), GBS, complete_reduction)
  GB.isGB = true
  GB.ord = ord_dest
  if isdefined(GB.gensBiPolyArray, :S)
    GB.gensBiPolyArray.S.isGB  = true
  end
  I.gb[ord_dest] = GB
  return GB
end

# Helper functions for groebner_basis_with_hilbert

function _extract_weights(T::MPolyDecRing)
  if !is_z_graded(T)
    error("Ring must be graded by the Integers.")
  end
  return [Int(first(gr_elem.coeff)) for gr_elem in T.d]
end

function _extend_mon_order(ordering::MonomialOrdering,
                           homogenized_ring::MPolyDecRing)

  nvars = ngens(ordering.R)
  m = canonical_matrix(ordering)
  m_hom = similar(m, nvars + 1, nvars + 1)
  m_hom[1, :] = ones(Int, nvars + 1)
  m_hom[2:end, 2:end] = m
  return matrix_ordering(homogenized_ring, m_hom)
end

function _mod_rand_prime(I::MPolyIdeal)
  p = 32771
  while true
    p = Hecke.next_prime(p)
    
    base_field = GF(p)
    ModP, _ = polynomial_ring(base_field, ngens(base_ring(I)); cached = false)
    I_mod_p_gens =
      try
        [map_coefficients(base_field, f; parent=ModP) for f in gens(I)]
      catch e
        # this precise error is thrown if the chosen prime p divides one
        # of the denominators of the coefficients of the generators of I.
        # In this case we simply choose the next prime and try again.
        if e == ErrorException("Unable to coerce") 
          continue
        else
          rethrow(e)
        end
      end
    return ideal(ModP, I_mod_p_gens)
  end
end

# check homogeneity w.r.t. some weights

function _is_homogeneous(f::MPolyRingElem, weights::Vector{Int})
  w = sum(weights .* first(exponents(f)))
  all(sum(weights .* e) == w for e in exponents(f))
end


# check homogeneity w.r.t. total degree
function _is_homogeneous(f::MPolyRingElem)
  leadexpv,tailexpvs = Iterators.peel(AbstractAlgebra.exponent_vectors(f))
  d = sum(leadexpv)
  for tailexpv in tailexpvs
    if d!=sum(tailexpv)
      return false
    end
  end
  return true
end

# compute weights such that F is a homogeneous system w.r.t. these weights
function _find_weights(F::Vector{P}) where {P <: MPolyRingElem}

  if all(_is_homogeneous, F)
    return ones(Int, ngens(parent(F[1])))
  end

  nrows = sum((length).(F)) - length(F)
  ncols = ngens(parent(first(F)))

  exp_diffs = permutedims(reduce(hcat, [e[i] - e[1] for e in
                                          (collect).((exponents).(F))
                                          for i in 2:length(e)]))
  K = kernel(matrix(QQ, nrows, ncols, exp_diffs); side = :right)
  isempty(K) && return zeros(Int, ncols)
  # Here we try to find a vector with strictly positive entries in K
  # this method to find such a vector is taken from
  # https://mathoverflow.net/questions/363181/intersection-of-a-vector-subspace-with-a-cone
  Pol = polyhedron(-K,  zeros(Int, ncols))
  !is_feasible(Pol) && return zeros(Int, ncols)
  pos_vec = [ZZ(0) for i in range(1,ncols)]
  for i in 1:ncols
    ei = [j == i ? one(QQ) : zero(QQ) for j in 1:ncols]
    obj_func = ei * K
    L = linear_program(Pol, obj_func)
    m, v = solve_lp(L)
    if isnothing(v)
      Pol_new = intersect(Pol, polyhedron(ei*K, [1]))
      L = linear_program(Pol_new, obj_func)
      v = optimal_vertex(L)
    end
    pos_vec += ((v.p)*transpose(K))[1,:]
  end
  ret = (Int).(lcm((denominator).(pos_vec)) .* pos_vec)
  ret = (x -> div(x, gcd(ret))).(ret) 
  # assure that the weights fit in Int32 for singular
  return all(ret .< 2^32) ? ret : zeros(Int,ncols)
end
