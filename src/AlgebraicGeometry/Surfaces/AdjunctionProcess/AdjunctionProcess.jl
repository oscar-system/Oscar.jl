
##############################################################################
#
# adjunction process for surfaces
#
##############################################################################

#function _random_vector(m::Int, I::MPolyIdeal, bound::Int = 30000)
#   @assert m>0
#   R = base_ring(I)  
#   SM = Singular.LibRandom.randommat(m, 1, Oscar.singular_generators(I), bound)
#   return [R(SM[i,1]) for i=1:nrows(SM)]
#end

#function _random_matrix(m::Int, n::Int, I::MPolyIdeal, bound::Int = 30000)
#   @assert n>0
#   @assert m>0
#   R = base_ring(I)  
#   SM = Singular.LibRandom.randommat(m, n, Oscar.singular_generators(I), bound)
#   return Matrix(matrix(R, nrows(SM), ncols(SM), [R(SM[i,j]) for i=1:nrows(SM) for j=1:ncols(SM)]))
#end

function _random_vector(R::MPolyDecRing, m::Int, d::Int)
   @assert m>0
   rd =  _random_matrix(R, m, 1, d)
   return [R(rd[i,1]) for i=1:nrows(rd)]
end

function _random_matrix(R::MPolyDecRing, m::Int, n::Int, d::Int)
   @assert n>0
   @assert m>0
   F = graded_free_module(R, m)
   G = graded_free_module(R, n)
   G = twist(G, d)
   H, h = hom(F, G)
   phi = rand_homogeneous(H, 0)
   return matrix(h(phi))
end

function _arithmetic_genus(I::MPolyIdeal)
   S = base_ring(I)
    #req is_standard_graded(S) "The base ring must be standard graded."
   A, _ = quo(S, I)
   H = hilbert_polynomial(A)
   return (-1)^(krull_dim(A)-1)*(ZZ(coeff(H, 0)) - 1)
end

@doc raw"""
    sectional_genus(X::AbsProjectiveVariety)

Given a subvariety `X` of some $\mathbb P^n$, return the arithmetic genus of the intersection of `X`
with a general linear subspace of $\mathbb P^n$ of dimension $c+1$.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> degree(X)
6

julia> sectional_genus(X)
3
```
"""
function sectional_genus(X::AbsProjectiveVariety)
   I = defining_ideal(X)
   S = base_ring(I)
   nl = length(gens(S))-codim(I)-2
   rd = _random_vector(S, nl, 1)
   J = ideal(S, rd)+I
   dd = dim(J)
   while !(dd == 2)
      rd = _random_vector(S, nl, 1)
      J = ideal(S, rd)+I
      dd = dim(J)
   end
   return _arithmetic_genus(J)
end

@doc raw"""
    is_linearly_normal(X::AbsProjectiveVariety)

Return `true` if `X` is linearly normal, and `false` otherwise.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> is_linearly_normal(X)
true
```
"""
function is_linearly_normal(X::AbsProjectiveVariety)
   I = defining_ideal(X)
   M = ideal_as_module(I)
   n = length(gens(base_ring(I)))
   tbl = sheaf_cohomology(M, -2, 1)
   return tbl[1, 1] == 0
end

function _adjoint_matrix(D::AbstractAlgebra.Generic.MatSpaceElem)
     #@req !any(x->is_zero(x) ? true : degree(Int, x) > 1, Oscar._vec(D)) "The code assumes that the input matrix has linear entries."
     R = base_ring(D)
     K = coefficient_ring(R)
     r = ncols(D)
     P, _ = graded_polynomial_ring(K, :z => 1:r; cached=false)
     RP, _ = graded_polynomial_ring(K, vcat(symbols(R), symbols(P)); cached=false)
     embRRP = hom(R, RP, gens(RP)[1:ngens(R)])    
     projRPP = hom(RP, P, vcat([zero(P) for i = 1:ngens(R)], gens(P)))
     M = map(embRRP, D)*matrix(gens(RP)[(ngens(R)+1):ngens(RP)])
     DM = [map_entries(x->derivative(x,i), M) for i = 1:ngens(R)]
     MM = transpose(hcat(DM...));
     return transpose(map(projRPP, MM))
end

function _matrix_has_no_quadratic_or_higher_entries(D::AbstractAlgebra.Generic.MatSpaceElem)
   return !any(x -> !is_zero(x) && total_degree(x) > 1, Oscar._vec(D))
end

function _is_zero_free_module_for_adjunction(M)
   try
      return rank(M) == 0
   catch
   end
   try
      return ngens(M) == 0
   catch
   end
   try
      return is_zero(M)
   catch
   end
   try
      return iszero(M)
   catch
   end
   return false
end

function _resolution_cached_zero_after_codimension(C, c::Int)
   try
      cf = chain_factory(C)
      if hasfield(typeof(cf), :map_cache) && length(cf.map_cache) >= c + 1
         return _is_zero_free_module_for_adjunction(domain(cf.map_cache[c + 1]))
      end
   catch
   end

   try
      cf = chain_factory(original_complex(C))
      if hasfield(typeof(cf), :map_cache) && length(cf.map_cache) >= c + 1
         return _is_zero_free_module_for_adjunction(domain(cf.map_cache[c + 1]))
      end
   catch
   end
   return false
end

function _canonical_presentation_matrix_via_module(X::AbsProjectiveVariety)
   Omega = canonical_bundle(X)
   FOmega = free_resolution(Omega, length = 1, algorithm = :mres)
   return matrix(map(FOmega, 1))
end

function _zero_presentation_matrix_for_adjunction(Pn::MPolyDecRing, ncols::Int = 1)
   return matrix(Pn, 0, ncols, elem_type(Pn)[])
end

function _minimal_resolution_for_adjunction(F)
   hasmethod(minimize, Tuple{typeof(F)}) && return minimize(F)
   return simplify(F)
end

function _canonical_presentation_matrix_direct(X::AbsProjectiveVariety)
   Pn = ambient_coordinate_ring(X)
   c = Int(codim(X))
   c == 0 && return _zero_presentation_matrix_for_adjunction(Pn)

   A = homogeneous_coordinate_ring(X)
   FA, _ = free_resolution(SimpleFreeResolution, A)
   C_min = _minimal_resolution_for_adjunction(FA)
   D = transpose(matrix(map(C_min, c)))
   _resolution_cached_zero_after_codimension(C_min, c) || return nothing
   return D
end

function _matrix_from_column_vectors_for_adjunction(K, cols::Vector)
   isempty(cols) && return matrix(K, 0, 0, elem_type(K)[])
   m = length(cols[1])
   return matrix(K, m, length(cols), [cols[j][i] for i in 1:m for j in 1:length(cols)])
end

function _coefficient_of_monomial_for_adjunction(f, m)
   e = exponent_vector(m, 1)
   for (c, ee) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
      ee == e && return c
   end
   return zero(coefficient_ring(parent(f)))
end

function _coordinates_in_quotient_degree_for_adjunction(f, projection_to_quotient, quotient_monomials::Vector)
   nf = lift(simplify(projection_to_quotient(f)))
   return [_coefficient_of_monomial_for_adjunction(nf, m) for m in quotient_monomials]
end

function _homogeneous_generators_for_adjunction(I::MPolyIdeal)
   result = elem_type(base_ring(I))[]
   for f in gens(I)
      is_zero(f) && continue
      if is_homogeneous(f)
         push!(result, f)
      else
         try
            append!(result, collect(values(homogeneous_components(f))))
         catch
            push!(result, f)
         end
      end
   end
   return result
end

function _degree_piece_candidates_of_ideal_for_adjunction(I::MPolyIdeal, d::Int)
   R = base_ring(I)
   d < 0 && return elem_type(R)[]

   result = elem_type(R)[]
   for f in _homogeneous_generators_for_adjunction(I)
      df = total_degree(f)
      df <= d || continue
      for m in monomial_basis(R, d - df)
         push!(result, m*f)
      end
   end
   return result
end

function _basis_of_ideal_quotient_degree_piece_for_adjunction(L::MPolyIdeal, J::MPolyIdeal, d::Int)
   R = base_ring(L)
   K = coefficient_ring(R)
   A, pi = quo(R, J)
   quotient_monomials = monomial_basis(A, d)

   basis_polys = elem_type(R)[]
   basis_vectors = Vector{Vector{elem_type(K)}}()
   current_rank = 0

   for f in _degree_piece_candidates_of_ideal_for_adjunction(L, d)
      nf = lift(simplify(pi(f)))
      v = [_coefficient_of_monomial_for_adjunction(nf, m) for m in quotient_monomials]
      any(x -> !is_zero(x), v) || continue

      test_vectors = copy(basis_vectors)
      push!(test_vectors, v)
      new_rank = rank(_matrix_from_column_vectors_for_adjunction(K, test_vectors))
      if new_rank > current_rank
         push!(basis_polys, nf)
         push!(basis_vectors, v)
         current_rank = new_rank
      end
   end

   return basis_polys
end

function _linear_relation_matrix_from_linkage_degree_piece_for_adjunction(R::MPolyDecRing, J::MPolyIdeal, basis_polys::Vector, d::Int)
   r = length(basis_polys)
   r == 0 && return _zero_presentation_matrix_for_adjunction(R, 0)

   K = coefficient_ring(R)
   A, pi = quo(R, J)
   target_monomials = monomial_basis(A, d + 1)
   vars = gens(R)
   nvars = length(vars)

   columns = Vector{Vector{elem_type(K)}}()
   for j in 1:r
      for i in 1:nvars
         push!(columns, _coordinates_in_quotient_degree_for_adjunction(vars[i]*basis_polys[j], pi, target_monomials))
      end
   end

   ker = kernel(_matrix_from_column_vectors_for_adjunction(K, columns), side = :right)
   nrels = ncols(ker)
   nrels == 0 && return _zero_presentation_matrix_for_adjunction(R, r)

   entries = elem_type(R)[]
   for a in 1:nrels
      for j in 1:r
         lin = zero(R)
         for i in 1:nvars
            lin += R(ker[(j - 1)*nvars + i, a])*vars[i]
         end
         push!(entries, lin)
      end
   end

   return matrix(R, nrels, r, entries)
end

function _try_complete_intersection_pair_for_adjunction(I::MPolyIdeal; random_trials::Int = 20)
   R = base_ring(I)
   G = sort(_homogeneous_generators_for_adjunction(I), by = total_degree)
   length(G) >= 2 || return nothing

   for i in 1:(length(G)-1)
      for j in (i+1):length(G)
         f = G[i]
         g = G[j]
         is_zero(f) && continue
         is_zero(g) && continue
         codim(ideal(R, [f, g])) == 2 && return (f, g)
      end
   end

   degrees = sort(unique(total_degree.(G)))
   for d in degrees
      Gd = [f for f in G if total_degree(f) == d]
      length(Gd) >= 2 || continue

      for _ in 1:random_trials
         C = _random_matrix(R, 2, length(Gd), 0)
         f = sum(C[1, j]*Gd[j] for j in 1:length(Gd); init = zero(R))
         g = sum(C[2, j]*Gd[j] for j in 1:length(Gd); init = zero(R))
         is_zero(f) && continue
         is_zero(g) && continue
         codim(ideal(R, [f, g])) == 2 && return (f, g)
      end
   end

   return nothing
end

function _canonical_presentation_matrix_via_ci_linkage_linear_strand(
    X::AbsProjectiveVariety;
    random_trials::Int = 20
  )
   R = ambient_coordinate_ring(X)
   is_standard_graded(R) || return nothing
   ngens(R) == 5 || return nothing
   Int(codim(X)) == 2 || return nothing

   I = defining_ideal(X)
   fg = _try_complete_intersection_pair_for_adjunction(I; random_trials)
   fg === nothing && return nothing

   f, g = fg
   J = ideal(R, [f, g])
   L = quotient(J, I)
   d = total_degree(f) + total_degree(g) - 4
   d < 0 && return _zero_presentation_matrix_for_adjunction(R, 0)

   basis_polys = _basis_of_ideal_quotient_degree_piece_for_adjunction(L, J, d)
   return _linear_relation_matrix_from_linkage_degree_piece_for_adjunction(R, J, basis_polys, d)
end

function _canonical_presentation_matrix_for_adjunction(X::AbsProjectiveVariety; algorithm::Symbol = :fast)
   if algorithm == :direct
      D = _canonical_presentation_matrix_direct(X)
      D !== nothing && return D
      error("The direct canonical presentation is not available for this embedding; use canonical_algorithm = :fast, :linkage_linear_strand, or :module.")
   elseif algorithm == :linkage_linear_strand
      D = _canonical_presentation_matrix_via_ci_linkage_linear_strand(X)
      D !== nothing && return D
      error("The complete-intersection linkage linear-strand presentation is only available for suitable codimension-two surfaces in projective 4-space.")
   elseif algorithm == :fast
      D = _canonical_presentation_matrix_direct(X)
      D !== nothing && return D
      D = _canonical_presentation_matrix_via_ci_linkage_linear_strand(X)
      D !== nothing && return D
      return _canonical_presentation_matrix_via_module(X)
   elseif algorithm == :module
      return _canonical_presentation_matrix_via_module(X)
   end
   error("Unsupported canonical presentation algorithm $(algorithm). Use :fast, :direct, :linkage_linear_strand, or :module.")
end

@doc raw"""
    canonical_bundle(X::AbsProjectiveVariety)

Given a smooth projective variety `X`, return a module whose sheafification is the canonical bundle of `X`.

!!! note
    The function does not check smoothness. If you are uncertain, enter `is_smooth(X)` first.

# Examples

```jldoctest
julia> R, x = graded_polynomial_ring(QQ, :x => (1:6));

julia> I = ideal(R, [x[1]*x[6] - x[2]*x[5] + x[3]*x[4]]);

julia> GRASSMANNIAN = variety(I);

julia> Omega = canonical_bundle(GRASSMANNIAN)
Graded subquotient of graded submodule of R^1 with 1 generator
  1: e[1]
by graded submodule of R^1 with 1 generator
  1: (x[1]*x[6] - x[2]*x[5] + x[3]*x[4])*e[1]

julia> degrees_of_generators(Omega)
1-element Vector{FinGenAbGroupElem}:
 [4]
```

```
julia> R, (x, y, z) = graded_polynomial_ring(QQ,[:x, :y, :z]);

julia> I = ideal(R, [y^2*z + x*y*z - x^3 - x*z^2 - z^3]);

julia> ELLCurve = variety(I);

julia> Omega = canonical_bundle(ELLCurve)
Graded submodule of R^1 with 1 generator
  1: e[1]
represented as subquotient with no relations

julia> degrees_of_generators(Omega)
1-element Vector{FinGenAbGroupElem}:
 [0]
```

```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> Omega = canonical_bundle(X);

julia> typeof(Omega)
SubquoModule{MPolyDecRingElem{fpFieldElem, fpMPolyRingElem}}
```
"""
function canonical_bundle(X::AbsProjectiveVariety)
  Pn = ambient_coordinate_ring(X)
  A = homogeneous_coordinate_ring(X)
  n = ngens(Pn)-1
  c = codim(X)
  FA, _ = free_resolution(SimpleFreeResolution, A)
  C_simp = simplify(FA)
  C_shift = shift(C_simp, c)
  OmegaPn = graded_free_module(Pn, [n+1])
  D = hom(C_shift, OmegaPn)
  D_simp = simplify(D)
  Z, inc = kernel(D_simp, 0)
  B, inc_B = boundary(D_simp, 0)
  return prune_with_map(SubquoModule(D_simp[0], ambient_representatives_generators(Z), ambient_representatives_generators(B)))[1]
end

@doc raw"""
    adjunction_process(X::AbsProjectiveVariety, steps::Int=0; canonical_algorithm::Symbol=:fast)

Given a smooth surface `X` and a non-negative integer `steps`, return data which describes the adjunction process for `X`: 
If `steps == 0`, carry out the complete process. Otherwise, carry out the indicated number of steps only.

The keyword argument `canonical_algorithm` controls how the presentation matrix of the canonical module is computed in each adjunction step. With the default `:fast`, OSCAR first tries to read this matrix directly from the codimension step of a minimal free resolution of the homogeneous coordinate ring. If this is certified, the matrix is returned even when it has quadratic or higher-degree entries; in that case the adjunction process stops, since the implemented adjoint matrix construction requires a linear presentation. If the direct extraction is not certified and the surface has codimension two in projective 4-space, OSCAR next tries a complete-intersection linkage computation of the linear strand of the canonical module. Only if both fast methods fail does `:fast` fall back to the generic route through `canonical_bundle`. The value `:module` yields the generic construction through `canonical_bundle` and a length-one free resolution of that module. The value `:direct` uses only the direct extraction and throws an error if the resolution does not certify that this extraction is valid. The value `:linkage_linear_strand` uses only the complete-intersection linkage linear-strand computation and throws an error if it is not applicable.

More precisely, if $X^{(0)} = X \rightarrow X^{(1)}\rightarrow \dots \rightarrow X^{(r)}$ is the sequence of successive adjunction maps and 
adjoint surfaces in the completed adjunction process, return a quadruple `L`, say, where: 

`L[1]` is a vector of tuples of numerical data: For each step $X^{(i)}\rightarrow X^{(i+1)}$, return the tuple $(n^{(i)}, d^{(i)}, \pi^{(i)}, s^{(i)}),$ 
where $n^{(i)}$ is the dimension of the ambient projective space of $X^{(i)}$, $d^{(i)}$ is the degree of $X^{(i)}$, $\pi^{(i)}$ is the sectional genus of $X^{(i)}$, 
and $s^{(i)}$ is the number of exceptional $(-1)$-lines on $X^{(i)}$ which are blown down to points in $ X^{(i+1)}$.

`L[2]` is a vector of adjoint matrices: For each step $X^{(i)}\rightarrow X^{(i+1)}$, return a presentation matrix of $S_X^{(i)}(1)$ considered 
as a module over $S_X^{(i+1)}$, where the $S_X^{(i)}$ are the homogeneous coordinate rings of the $X^{(i)}$. If `X` is rational, these matrices
can be used to compute a rational parametrization of `X`.

`L[3]` is a vector of zero-dimensional projective algebraic sets: For each step $X^{(i)}\rightarrow X^{(i+1)}$, return the union of points in $ X^{(i+1)}$
which are obtained by blowing down the exceptional $(-1)$-lines on $X^{(i)}$.
    
`L[4]` is a projective variety: Return the last adjoint surface $X^{(r)}$. 


!!! note
    The function does not check whether `X` is smooth. If you are uncertain, enter `is_smooth(X)` first.

!!! warning
    At current state, the adjunction process is only implemented for rational and Enriques surfaces which are linearly normal in the given embedding. The function does not check whether `X` is rational or an Enriques surface. In fact, at current state, `OSCAR` does not offer direct checks for this. Note, however, that the adjunction process will give an answer to this question a posteriori in cases where it terminates with a surface which is known to be rational or an Enriques surface.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators


julia> dim(X)
2

julia> codim(X)
2

julia> L = adjunction_process(X);

julia> L[1]
2-element Vector{NTuple{4, ZZRingElem}}:
 (4, 6, 3, 0)
 (2, 1, 0, 10)

julia> L[4]
Projective variety
  in projective 2-space over GF(31991) with coordinates [z[1], z[2], z[3]]
defined by ideal (0)

julia> L[3][1]
Projective algebraic set
  in projective 2-space over GF(31991) with coordinates [z[1], z[2], z[3]]
defined by ideal with 5 generators

julia> dim(L[3][1])
0

julia> degree(L[3][1])
10
```

```
julia> X = rational_d9_pi7();

julia> L = adjunction_process(X);

julia> L[1]
3-element Vector{NTuple{4, ZZRingElem}}:
 (4, 9, 7, 0)
 (6, 9, 4, 6)
 (3, 3, 1, 3)
```

!!! note
    Inspecting the  returned numerical data in the first example above, we see that the Bordiga surface is the blow-up of the projective plane in 10 points, embedded into projective 4-space by the linear system $H = 4L -\sum_{i=1}^{10} E_i$. Here, $L$ is the preimage of a line and the $E_i$ are the exceptional divisors. In the second example, we see from the output that the terminal object of the adjunction process is a Del Pezzo surface in projective 3-space, that is, the blow-up of the projective plane in 6 points. In sum, we see that `X` is the blow-up of the projective plane in 15 points, embedded into projective 4-space by the linear system $H = 9L - \sum_{i=1}^{6} 3E_i - \sum_{i=7}^{9} 2E_i - \sum_{i=10}^{15} E_i$. 
"""
function adjunction_process(X::AbsProjectiveVariety, steps::Int = 0; canonical_algorithm::Symbol = :fast)
   @assert steps >= 0
   Pn = ambient_coordinate_ring(X)
   I = defining_ideal(X)
   Y = X
   @req dim(I) == 3 "The given variety is not a surface."
   @req is_linearly_normal(X) "The given variety is not linearly normal."
   # TODO If X is not linear normal, embed it as a linearly normal variety first.
   numlist = [(ZZ(ngens(Pn)-1), degree(X), sectional_genus(X), zero(ZZ))]
   ptslist = ProjectiveAlgebraicSet[]
   adjlist = AbstractAlgebra.Generic.MatSpaceElem[] 
   D = _canonical_presentation_matrix_for_adjunction(X; algorithm = canonical_algorithm)
   count = 1
   while nrows(D) > 2 && (steps == 0 || count <= steps)
      if _matrix_has_no_quadratic_or_higher_entries(D)
         adj = _adjoint_matrix(D)
      else
        return (numlist, adjlist, ptslist, variety(I, check = false, is_radical = true))
      end
      rd = _random_matrix(Pn, 3, ncols(adj), 0)
      dd = dim(ideal(Pn, rd*gens(Pn))+I) 
      while !(dd == 0)
         rd = _random_matrix(Pn, 3, ncols(adj), 0)
         dd = dim(ideal(Pn, rd*gens(Pn))+I)
      end

      I = annihilator(cokernel(adj))
      Pn = base_ring(I)
    
      Ipts = ideal(Pn, [zero(Pn)])    
      RD =  map_entries(constant_coefficient, rd)   
      for i = 1:3
          KKi = map_entries(Pn, kernel(matrix(RD[i,:])))
          AAi = annihilator(cokernel(adj*transpose(KKi)))
	  Ipts = Ipts+AAi
      end
      Ipts = saturation(Ipts, ideal(Pn, gens(Pn)))
      pts = algebraic_set(Ipts,  is_radical = true, check = false) 
      if dim(pts) == 0
         l = degree(pts)
      else
         l = zero(ZZ)
      end
      Y = variety(I, check = false, is_radical = true)
      dY = degree(Y)
      piY = sectional_genus(Y)
      dummy = (ZZ(ngens(Pn)-1), dY, piY)
      if l==0 && dummy == numlist[count][1:3]   # Enriques surface
	  return (numlist, adjlist, ptslist, Y)
      else
        push!(numlist, (ZZ(ngens(Pn)-1), dY, piY, l))
        push!(adjlist, adj)
        push!(ptslist, pts)
        D = _canonical_presentation_matrix_for_adjunction(Y; algorithm = canonical_algorithm)
      end
      count = count+1
   end
   return (numlist, adjlist, ptslist, Y)
end
