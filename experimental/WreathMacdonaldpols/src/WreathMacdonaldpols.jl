# I would like to credit Dario Mathiä who produced an initial version of the following code

export wreath_macdonald_polynomials, wreath_macdonald_polynomial

# Tools

# Computes the b-invariant
function b1(lambda::Vector{Int})
  return sum((i-1)*l_i for (i, l_i) in enumerate(lambda); init=0)
end

# Computes the b-invariant of a multipartition
function b_inv(lbb::Multipartition)
  return b1(map(sum,lbb)) + length(lbb)*sum(lambda -> b1(data(lambda)),lbb)
end

function beta_number_to_partition(beta::Vector{Int})
  lb=Int[]
  for j in 1:length(beta)
    nbholes=beta[j]-beta[1]-(j-1)
    if nbholes >= 1
      pushfirst!(lb, nbholes)
    end
  end
  return partition(lb)
end

# Computes the core associated to the coroot
#The charge of the core is equal to the coroot
#in the canonical basis (e_i) of Z^I.
function core(coroot::Vector{Int},r::Int)
  beta=Int[]
  m=minimum(coroot)
  M=maximum(coroot)
  beta=[k*r+(j-1) for k in m:M for j in 1:r if k <= coroot[j]]
  return beta_number_to_partition(beta)
end

#Inspired by Sagemath's algorithm to obtain a partition from core and quotient
#Be careful: it is not the usual quotient but the one defined in [Gordon 2008]
function tau_om(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int})
   r=length(lbb)
   w0=perm([r-i for i in 0:(r-1)])  #cf Remark 3.5 Orr and Shimozono
   lbb_perm=multipartition(permuted(lbb.mp,wperm*w0))
   gamma=core(coroot,r)
   lg=length(gamma)
   k=r*maximum(length(lbb_perm[i]) for i in 1:r) + lg
   v=[[gamma[i]-(i-1) for i in 1:lg]; [-i for i in lg:(k-1)]]
   w=[[x for x in v if mod((x-i),r) == 0] for i in 1:r]
   new_w=Int[]
   for i in 1:r
     lw=length(w[i])
     lq=length(lbb_perm[i])
     append!(new_w, w[i][1:lq] + r*lbb_perm[i])
     append!(new_w, w[i][lq+1:lw])
   end
   sort!(new_w,rev=true)
   new_w=[new_w[i]+(i-1) for i in 1:length(new_w)]
   filter!(x-> x!=0,new_w)
   return partition(new_w)
end

#Reoders the CharTable from Chevie
function reorder(charTirr::Matrix{QQAbFieldElem{AbsSimpleNumFieldElem}}, Modules::Vector{Multipartition{Int}},  mps::Vector{Multipartition{Int}})
  new_ord=indexin(mps,Modules)
  return charTirr[new_ord,:]
end

#Tools from representation theory

# Computes the fake degree of a multipartition cf. Ste89 Thm 5.3 and Prop. 3.3.2 Haiman cdm
function fake_deg(lbb::Multipartition, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}}, var::AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}})
  r=length(lbb)
  n=sum(lbb)
  res=Q(1)
  for p in 1:n
    res=res*(1-var^(r*p))
  end
  for lambda in lbb
    for i in 1:length(lambda)
      for j in 1:lambda[i]
        hookfactor=(Q(1)-var^(r*(1+lambda[i]+count(>=(j),lambda)-j-i)))
        res=res//hookfactor
      end
    end
  end
  res=res*var^b_inv(lbb)
  return res
end

# computes C_Delta defined in PhD relation (1.45)
function C_Delta(r::Int, n::Int, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}}, var::AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}, mps::Vector{Multipartition{Int}})
  charTable=character_table_complex_reflection_group(r,1,n)
  charTirr=[charTable[i,j] for i in 1:nrows(charTable), j in 1:ncols(charTable)]
  Modules=[multipartition([lbb...]) for lbb in class_parameters(charTable)]
  l=length(Modules)
  charTirr=reorder(charTirr,Modules,mps)
  charTirrT=solve_init(matrix(Q,transpose(charTirr)))
  rows=zero_matrix(Q,l,0)
  for i in 1:l
    col=zero_matrix(Q,l,1)
    for j in 1:l
      v=matrix(Q,l,1,map(k->charTirr[i,k]*charTirr[j,k],1:l))
      x=solve(charTirrT,v,side=:right)
      col=col+fake_deg(mps[j],Q,var)*x
    end
    rows=hcat(rows,col)
  end
  return rows
end

# computes the character of the simple representations of the Cherednik algebra.
function C_L(r::Int, n::Int, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}}, var::AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}, mps::Vector{Multipartition{Int}})
  C_D=C_Delta(r,n,Q,var,mps)
  d = [var^b_inv(mp)//fake_deg(mp,Q,var) for mp in mps]
  diag=diagonal_matrix(Q,d)
  return diag * C_D
end

function bigger_ord(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int}, mps::Vector{Multipartition{Int}})
  r=length(lbb)
  n=sum(lbb)
  lbbquot=tau_om(lbb,wperm,coroot)
  res=map(x->x,Iterators.filter(x-> dominates(tau_om(x,wperm,coroot),lbbquot),mps))
  return res
end

function smaller_ord(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int}, mps::Vector{Multipartition{Int}})
  r=length(lbb)
  n=sum(lbb)
  lbbquot=tau_om(lbb,wperm,coroot)
  res=map(x->x,Iterators.filter(x-> dominates(lbbquot,tau_om(x,wperm,coroot)),mps))
  return res
end

@doc raw"""
    wreath_macdonald_polynomials(n::Int,
                                 r::Int,
                                 wperm::PermGroupElem,
                                 coroot::Vector{Int},
                                 parent::AbstractAlgebra.Generic.MPolyRing{QQAbFieldElem{AbsSimpleNumFieldElem}})

Given two integers n and r and an element of the affine Weyl group of type A (seen as
the semi-direct product of the Symmetric group with the coroot lattice), this function
returns the square matrix of coefficients of the wreath Macdonald polynomials associated with
all multipartition of size n and length r in the standard Schur basis indexed by the
multipartitions of the same size and length as lbb. Each row of this matrix is a wreath
Macdonald polynomial.
"""
function wreath_macdonald_polynomials(n::Int,
                                      r::Int,
                                      wperm::PermGroupElem,
                                      coroot::Vector{Int};
                                      parent::AbstractAlgebra.Generic.MPolyRing{QQAbFieldElem{AbsSimpleNumFieldElem}})

  Q = fraction_field(parent)
  q,t = gens(parent)

  mps=collect(multipartitions(n,r))
  l=length(mps)

  c_L=C_L(r,n,Q,t,mps)
  c_L_q=map_entries(f->numerator(f)(0,q)//denominator(f)(0,q),c_L)
  c_L_tinv=map_entries(f->numerator(f)(0,1//t)//denominator(f)(0,1//t),c_L)

  rows=[]
  for lbb in mps
    smallers=smaller_ord(lbb, wperm, coroot,mps) #tinv
    biggers=bigger_ord(lbb, wperm, coroot,mps) #q
    smaller_indices=sort!(indexin(smallers,mps))
    bigger_indices=sort!(indexin(biggers,mps))
    sub_smaller_tinv=vcat(map(i->c_L_tinv[i:i,:],smaller_indices)...)
    sub_bigger_q=vcat(map(i->c_L_q[i:i,:],bigger_indices)...)
    M=vcat(sub_smaller_tinv,sub_bigger_q)
    B=kernel(M, side=:left)
    rows=vcat(rows,B[1,1:length(smaller_indices)]*sub_smaller_tinv)
  end
  c_L_qt = matrix(Q,l,l,rows)

  triv=[partition([n])]
  for i in 2:r
    push!(triv,partition([]))
  end
  triv=multipartition(triv)
  index_triv=findfirst(x-> x==triv,mps)
  d=[1//c_L_qt[i,index_triv] for i in 1:l]
  diag=diagonal_matrix(Q,d)
  c_L_qt_H=diag*c_L_qt
  return c_L_qt_H
end

@doc raw"""

    wreath_macdonald_polynomial(lbb::Multipartition,
                                wperm::PermGroupElem,
                                coroot::Vector{Int},
                                parent::AbstractAlgebra.Generic.MPolyRing{QQAbFieldElem{AbsSimpleNumFieldElem}})

Given a multipartition lbb and an element of the affine Weyl group of type A (seen as
the semi-direct product of the Symmetric group with the coroot lattice), this function
returns the coefficients of the wreath Macdonald polynomial associated with lbb and
the affine Weyl group element in the standard Schur basis indexed by the multipartitions
 of the same size and length as lbb.
"""
function wreath_macdonald_polynomial(lbb::Multipartition,
                                     wperm::PermGroupElem,
                                     coroot::Vector{Int};
                                     parent::AbstractAlgebra.Generic.MPolyRing{QQAbFieldElem{AbsSimpleNumFieldElem}})

  Q = fraction_field(parent)
  q,t = gens(parent)

  r=length(lbb)
  n=sum(lbb)
  mps=collect(multipartitions(n,r))
  l=length(mps)

  c_L=C_L(r,n,Q,t,mps)
  c_L_q=map_entries(f->numerator(f)(0,q)//denominator(f)(0,q),c_L)
  c_L_tinv=map_entries(f->numerator(f)(0,1//t)//denominator(f)(0,1//t),c_L)

  smallers=smaller_ord(lbb, wperm, coroot,mps) #tinv
  biggers=bigger_ord(lbb, wperm, coroot,mps) #q
  smaller_indices=sort!(indexin(smallers,mps))
  bigger_indices=sort!(indexin(biggers,mps))
  sub_smaller_tinv=vcat(map(i->c_L_tinv[i:i,:],smaller_indices)...)
  sub_bigger_q=vcat(map(i->c_L_q[i:i,:],bigger_indices)...)
  M=vcat(sub_smaller_tinv,sub_bigger_q)
  B=kernel(M, side=:left)
  row=B[1,1:length(smaller_indices)]*sub_smaller_tinv
  c_L_qt = matrix(Q,1,l,row)

  triv=[partition([n])]
  for i in 2:r
    push!(triv,partition([]))
  end
  triv=multipartition(triv)
  index_triv=findfirst(x-> x==triv,mps)
  c_L_qt_H=1//c_L_qt[index_triv]*c_L_qt
  return c_L_qt_H
end
