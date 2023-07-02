@doc raw"""
    has_du_val_singularities(X::Scheme)

Return whether the given `X` has at most du Val (surface) singularities.

# Example:
```jldoctest
julia> R,(x,y,z,w) = QQ["x","y","z","w"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x, y, z, w])

julia> I = ideal(R,[w,x^2+y^3+z^4])
ideal(w, x^2 + y^3 + z^4)

julia> Rq,_ = quo(R,I)
(Quotient of multivariate polynomial ring by ideal with 2 generators, Map from
Multivariate polynomial ring in 4 variables over QQ to Rq defined by a julia-function with inverse)

julia> X = Spec(Rq)
Spec of Quotient of multivariate polynomial ring by ideal with 2 generators

julia> has_du_val_singularities(X)
true

```
"""
function has_du_val_singularities(X::AbsProjectiveScheme{<:Field,<:Any})
  return has_du_val_singularities(covered_scheme(X))
end

function has_du_val_singularities(X::AbsSpec{<:Field,<:Any})
  R = OO(X)
  I = modulus(R)
  
  J = image_ideal(singular_locus(X)[2])
  J = ideal(base_ring(I), lift.(gens(J))) + I
  dim(J) == 0 || return false                            ## non-isolated
  return is_du_val_singularity(X,J)
end

function has_du_val_singularities(X::AbsCoveredScheme{<:Field})
  C = (has_attribute(X0, :simplified_covering) ? simplified_covering(X0) : default_covering(X0))

  I = ideal_sheaf_of_singular_locus(X)
  decomp = minimal_associated_points(I)

  ## we do the double loop here to avoid unnecessary checks
  for J in decomp
    for U in C
      !isone(J(U)) || continue
      is_du_val_singularity(X(U),J(U)) || return false    ## testing the point in one chart suffices
      break                         
    end
  end

  return true
end

@doc raw"""
    is_du_val_singularity(X::AbsSpec, I::Ideal)

Return whether the given X has at most du Val (surface) singularities at the geometric points
specified by the ideal I.

**Note**: For the ideal I in a ring R, dim(R/I) = 0 is asserted

# Example:
```jldoctest
julia> R,(x,y,z,w) = QQ["x","y","z","w"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x, y, z, w])

julia> I = ideal(R,[w,x^2+y^3+z^4])
ideal(w, x^2 + y^3 + z^4)

julia> Rq,_ = quo(R,I)
(Quotient of multivariate polynomial ring by ideal with 2 generators, Map from
Multivariate polynomial ring in 4 variables over QQ to Rq defined by a julia-function with inverse)

julia> J = ideal(R,[x,y,z,w])
ideal(x, y, z, w)

julia> X = Spec(Rq)
Spec of Quotient of multivariate polynomial ring by ideal with 2 generators

julia> is_du_val_singularity(X,J)
true

```
"""
function is_du_val_singularity(X::AbsSpec{<:Field,<:Any},I::Ideal)
  OOX = OO(X)
  dim(X) == 2 || error("Scheme not of dimension 2")
  J = modulus(OOX)
  !isone(J) || error("Scheme is empty")
  !iszero(J) || return true                  ## X smooth

  R = base_ring(I)
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")
  base_ring(OOX) === R || error("base rings need to be identical")

  dim(I) == 0 || error("second argument does not describe a 0-dimensional scheme")
  if get_attribute(I,:is_absolutely_prime, false)
    return _check_duval_at_point(J,I)[1]
  end

  decomp = absolute_primary_decomposition(I)

  for (_,_,I2,mult) in decomp
    set_attribute!(I2,:is_absolutely_prime,true)
    ## pass to algebraic extension
    r_changed = base_ring(I2)
    kk = coefficient_ring(r_changed)
    J_changed = ideal(r_changed,  [change_coefficient_ring(kk,a, parent = r_changed) for a=gens(J)])
    is_du_val_singularity(Spec(quo(r_changed,J_changed)[1]),I2) || return false
  end
  
  return true
end

@doc raw"""
    decide_du_val_singularity(X::AbsSpec, I::Ideal)

Return a vector of tuples T with the following data:
- T[1]::Bool answers 'X' has at most du Val (surface) singularities at the geometric points
specified by the ideal 'I'.
- T[2]::Ideal I_P associated prime of I (possibly over a suitable field extension)
  describing some geometrically irreducible point
- T[3]::Tuple type of the singularity at P
- T[4]::Int number of conjugate points

If X has a least one singularity which is not du Val, the returned vector contains a
single tuple T, with the following values:
- T[1] = false
- T[2] = point at which some non-du-Val singularity is present
- T[3] = empty tuple
- T[4] = 1

**Note**: For the ideal I in a ring R, dim(R/I) = 0 is asserted

# Example:
```jldoctest
julia> R,(x,y,z,w) = QQ["x","y","z","w"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x, y, z, w])

julia> I = ideal(R,[w,x^2+y^3+z^4])
ideal(w, x^2 + y^3 + z^4)

julia> Rq,_ = quo(R,I)
(Quotient of multivariate polynomial ring by ideal with 2 generators, Map from
Multivariate polynomial ring in 4 variables over QQ to Rq defined by a julia-function with inverse)

julia> J = ideal(R,[x,y,z,w])
ideal(x, y, z, w)

julia> decide_du_val_singularity(X,J)
1-element Vector{Tuple{Bool, MPolyIdeal{QQMPolyRingElem}, Tuple{Symbol, Int64}, Int64}}:
 (1, ideal(x, y, z, w), (:E, 6), 1)

""" 
function decide_du_val_singularity(X::AbsSpec{<:Field,<:Any},I::MPolyIdeal)
  OOX = OO(X)
  dim(X) == 2 || error("Scheme not of dimension 2")
  J = modulus(OOX)
  !isone(J) || error("Scheme is empty")
  !iszero(J) || return true                  ## X smooth

  R = base_ring(I)
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")
  base_ring(OOX) === R || error("base rings need to be identical")

  dim(I) == 0 || error("second argument is not a point")
  if get_attribute(I,:is_absolutely_prime, false)
    li = _check_duval_at_point(J,I)
    return [(li[1], I, li[2],1)]
  end

  result_vector = []
  decomp = absolute_primary_decomposition(I)
  for (_,_,I2,mult) in decomp
    set_attribute!(I,:is_absolutely_prime,true)
    ## pass to algebraic extension
    r_changed = base_ring(I2)
    kk = coefficient_ring(r_changed)
    J_changed = ideal(r_changed,  [change_coefficient_ring(kk,a, parent = r_changed) for a=gens(J)])
    tempvec = decide_du_val_singularity(Spec(quo(r_changed,J_changed)[1]),I2)
    for x in tempvec
      x[1] || return [x]
      x[4] = mult
      push!(result_vector,x)
    end
  end

  return result_vector
end

@doc raw"""
    _check_du_val_at_point(IX:Ideal,Ipt::Ideal))

Returns a tuple T with the following data:
- T[1]::Bool has V(IX) at most a du Val singularity at V(Ipt)
- T[2]::Tuple Type of du Val singularity at V(Ipt)

**Note**: Assumes Ipt to be absolutely irreducible.

**Note**: Do not call directly, only via the higher level functions is_du_Val_singularity and decide_du_Val_singularity.

"""
function _check_duval_at_point(IX::Ideal,Ipt::Ideal)
  R = base_ring(IX)
  R == base_ring(Ipt) || error("basering mismatch")
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")
 
  JM = jacobi_matrix(gens(IX))

  smooth = (:A,0)

  ## localize at point
  a = Oscar.rational_point_coordinates(Ipt)
  U = complement_of_point_ideal(R,a)
  RL, loc_map = localization(R,U)
  IX_loc = loc_map(IX)
  JM_loc =  map(x ->loc_map(x), JM[:,:])

  if !all(iszero(a))
    F_loc = free_module(RL,ngens(IX))
    Jm = sub(F_loc,JM_loc)
    Jm = Jm + loc_map(IX)*F_loc
    Jm_shifted = Oscar. shifted_module(Jm)
    F_shifted = ambient_free_module(Jm_shifted)
  else
    F = free_module(R,ngens(IX))
    Jm = sub(F,JM)[1]
    Jm_shifted = Jm + (IX * F)[1]
    F_shifted = F
  end

  o = negdegrevlex(R)*lex(F_shifted)
  F1 = leading_module(Jm_shifted,o)
  F1quo = quo(F_shifted, F1)[1]

  constant_mons = vector_space_dimension(F1quo,0) 
  constant_mons < 2 || return (false, typeof(smooth))                   ## not a hypersurface
  constant_mons > 0 || return (true, smooth)                            ## no singularity
  
  tau = vector_space_dimension(F1quo)

  corank = vector_space_dimension(F1quo,1)
  corank < 3 || return (false,typeof(smooth))                           ## at least T_3,3,3 not duVal
  corank > 1 || return (true,(:A,tau))                                  ## A_series

  # we now already know essentially a hypersurface of corank 2, count degrees of freedom cubicpart
  cubiccount = vector_space_dimension(F1quo,2)
  cubiccount < 3 || return  (false, typof(smooth))                      ## at least X_9
  cubiccount > 1 || return  (true, (:D,tau))                            ## D_series

  # it is definitely in the E/J series
  tau < 9 || return(false, typeof(smooth))                              ## at least J_10
  return (true, (:E,tau))                                               ## E_6, E_7, E_8
end

@doc raw"""
    vector_space_dimension(M::SubquoModule, d::Int)

Let R be a MPolyAnyRing over a field kk and let M be a subquotient module over R.
Then the command returns the dimension of the kk-vectorspace corresponding to the
degree d slice of M, where the degree of each variable of R is counted as one and
the one of each generator of the ambient free module of M as zero.

    vector_space_dimension(M::SubquoModule)

If M happens to be finite-dimensional as a kk-vectorspace, this returns its dimension.

# Examples:
```jldoctest
julia> R,(x,y,z,w) = QQ["x","y","z","w"];

julia> F = free_module(R,2);

julia> M,_ = quo(F,[1*gen(F,1),x^2*gen(F,2),y^3*gen(F,2),z*gen(F,2),w*gen(F,2)]);

julia> vector_space_dimension(M,1)
2

julia> vector_space_dimension(M,2)
2

julia> vector_space_dimension(M,3)
1

julia> vector_space_dimension(M)
6

```
"""
function vector_space_dimension(M::SubquoModule)
  
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  has_monomials_on_all_axes(LM) || error("not a finite dimensional vector space")
  
  d = 0
  sum_dim = 0
  tempdim = vector_space_dimension(M,0)

  while tempdim > 0
    sum_dim = sum_dim + tempdim
    d = d+1
    tempdim = vector_space_dimension(M,d)
  end
 
  return sum_dim
end

function vector_space_dimension(M::SubquoModule,d::Int64)
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  return length([x*e for x in Oscar.all_monomials(R, d) for e in gens(F) if !(x*e in LM)])
end
  
function vector_space_dimension(M::SubquoModule{T}
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  M_shift,_,_ = shifted_module(M)
  o = negdegrevlex(base_ring(M_shift))*lex(ambient_free_module(M_shift)
  LM = leading_module(M_shift,o)
  return vector_space_dimension(LM)
end

function vector_space_dimension(M::SubquoModule{T},d::Int64
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  M_shift,_,_ = shifted_module(M)
  o = negdegrevlex(base_ring(M_shift))*lex(ambient_free_module(M_shift)
  LM = leading_module(M_shift,o)
  return vector_space_dimension(LM,d)
end

function vector_space_dimension(M::SubquoModule{T}
  ) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

function vector_space_dimension(M::SubquoModule{T},d::Int64
  ) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

@doc raw"""
    vector_space_basis(M::SubquoModule, d::Int)

Let R be a MPolyAnyRing over a field kk and let M be a subquotient module over R.
Then the command returns a monomial basis of the kk-vectorspace corresponding to the
degree d slice of M, where the degree of each generator of R is counted as one and
the one of each generator of the ambient free module of M as zero.

    vector_space_basis(M::SubquoModule)

If M happens to be finite-dimensional as a kk-vectorspace, this returns a monomial basis of it.

# Examples:
```jldoctest
julia> R,(x,y,z,w) = QQ["x","y","z","w"];

julia> F = free_module(R,2);

julia> M,_ = quo(F,[1*gen(F,1),x^2*gen(F,2),y^3*gen(F,2),z*gen(F,2),w*gen(F,2)]);

julia> vector_space_basis(M,2)
2-element Vector{FreeModElem{QQMPolyRingElem}}:
 x*y*e[2]
 y^2*e[2]

julia> vector_space_basis(M,0)
1-element Vector{FreeModElem{QQMPolyRingElem}}:
 e[2]

julia> vector_space_basis(M)
6-element Vector{Any}:
 e[2]
 x*e[2]
 y*e[2]
 x*y*e[2]
 y^2*e[2]
 x*y^2*e[2]

```
"""
function vector_space_basis(M::SubquoModule)
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  has_monomials_on_all_axes(LM) || error("not a finite dimensional vector space")
  
  d = 0
  all_mons=[]
  temp_mons = vector_space_basis(M,0)

  while length(temp_mons) > 0
    append!(all_mons,temp_mons)
    d = d+1
    temp_mons=vector_space_basis(M,d)
  end

  return all_mons
end

function vector_space_basis(M::SubquoModule,d::Int64)
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  return [x*e for x in Oscar.all_monomials(R, d) for e in gens(F) if !(x*e in LM)]
end

function vector_space_basis(M::SubquoModule{T}
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  M_shift,_,_ = shifted_module(M)
  return vector_space_basis(M_shift)
end

function vector_space_basis(M::SubquoModule{T},d::Int64
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  M_shift,_,_ = shifted_module(M)
  return vector_space_basis(M_shift,d)
end

function vector_space_basis(M::SubquoModule{T}
  ) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

function vector_space_basis(M::SubquoModule{T},d::Int64
  ) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

@doc raw"""
    has_monomials_on_all_axes(M::SubquoModule)

Internal function to test whether M is finite-dimensional vector space. Do not use directly
"""
function has_monomials_on_all_axes(M::SubquoModule)
  R = base_ring(M)

  length(rels(M)) == 0 || error("not implemented for quotients")
  
  ambient_rank = ngens(ambient_free_module(M))
  genlist = ambient_representatives_generators(M)
  explist = Tuple{Vector{Int64}, Int64, Int}[]
  for x in genlist
    tempexp = leading_exponent(x)
    tempdeg = sum(tempexp[1])
    push!(explist,(tempexp[1],tempexp[2],tempdeg))
  end
  all(findfirst(x -> (x[1][i] == x[3] && x[2]==j),explist) != nothing for i in 1:ngens(R) for j in 1:ambient_rank) || return false
  return true
end

