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
    All weights must be positive. If no weight vector is entered by the user, all weights 
    are set to 1. An error is thrown if the generators of `I` are not homogeneous with 
    respect to the corresponding (weighted) degree.  

!!! note
    If $R$ denotes the parent ring of $I$, and $p, q\in\mathbb Z[t]$ are polynomials
    such that $p/q$ represents the Hilbert series of $R/I$ as a rational function with 
    denominator $q = (1-t^{w_1})\cdots (1-t^{w_n}),$ where $n$ is the number of variables 
    of $R$, and $w_1, \dots, w_n$ are the assigned weights, then `hilbert_numerator` is 
    meant to be $p$. If this numerator is not entered by the user, it will be computed 
    internally.

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
  
  all(f -> _is_homogeneous(f, weights), gens(I)) || error("I must be given by generators homogeneous with respect to the given weights.")
  isa(coefficient_ring(I), AbstractAlgebra.Field) || error("The underlying coefficient ring of I must be a field.")
  ordering = destination_ordering
  is_global(ordering) || error("Destination ordering must be global.")
  haskey(I.gb, ordering) && return I.gb[ordering]
  if isnothing(hilbert_numerator)
    if isempty(I.gb)
      J = iszero(characteristic(base_ring(I))) ? _mod_rand_prime(I) : I
      G = standard_basis(J, ordering=wdegrevlex(base_ring(J), weights))
    else
      G = standard_basis(I)
    end

    if characteristic(base_ring(I)) > 0 && ordering == wdegrevlex(base_ring(I), weights)
      return G
    end
    h = Singular.hilbert_series(singular_generators(G, G.ord), weights)

  else
    # Quoting from the documentation of Singular.hilbert_series:
    # The coefficient vector is returned as a `Vector{Int32}`, and the last element is not actually part of the coefficients of Q(t).
    # what?
    h = (Int32).([coeff(hilbert_numerator, i) for i in 0:degree(hilbert_numerator)+1])
  end

  singular_I_gens = singular_generators(I.gens, ordering)
  singular_ring = base_ring(singular_I_gens)
  J = Singular.Ideal(singular_ring, gens(singular_I_gens)...)
  i  = Singular.std_hilbert(J, h, (Int32).(weights),
                            complete_reduction = complete_reduction)
  GB = IdealGens(I.gens.Ox, i, complete_reduction)
  GB.isGB = true
  GB.ord = ordering
  if isdefined(GB, :S)
    GB.S.isGB  = true
  end
  I.gb[destination_ordering] = GB
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
  pos_vec = zeros(Int, ncols)
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
    pos_vec += K*(v.p)
  end
  ret = (Int).(lcm((denominator).(pos_vec)) .* pos_vec)
  ret = (x -> div(x, gcd(ret))).(ret) 
  # assure that the weights fit in Int32 for singular
  return all(ret .< 2^32) ? ret : zeros(Int,ncols)
end
