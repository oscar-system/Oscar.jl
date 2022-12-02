
######################################### ZERODIMENSIONAL PRIMARY DECOMPOSITION VIA GTZ #################################################################
# This is an implementation of the GTZ Algorithm, primarily based on Chapter 4 of "A SINGULAR introduction to Commutative Algebra" by Greuel and Pfister.
#
# Main algorithms:
# zerodecomp   
# primaryTest
#########################################################################################################################################################


######################################### zerodecomp ####################################################################################################

#Implementation of Algorithm 4.2.7
@doc Markdown.doc"""
    zerodecomp(I::Singular.sideal, outputReduced::Bool = false)

Compute the primary decomposition of a zero-dimensional Ideal $I$ over a basefield of characteristig 0 that has a global ordering,
via GTZ. Returns the primary decomposition with associated primes in a list. If $outputReduced$ is set to true, calculates a reduced Gröbner
before returning the primary decomposition.
"""
#REQUIRES:  INPUT sideal I, Optional: Whether the output should be put in completely reduced Gröbner bases. Default is to not reduce.
#           I is zero-dimensional (Krull dimension)  
#           I has a field as coefficient ring
#           basering is not a quotient ring
#           basering has a global ordering
#           basering has characteristic 0
#OUTPUT:    Primary decomposition of I with associated primes in a list 

function zerodecomp(I::Singular.sideal, outputReduced::Bool = false, usefglm::Bool = false)
  if ordering(I.base_ring) != :lex
     Icopy = changeOrderOfBasering(I, :lex);
  else 
     Icopy = Singular.deepcopy(I)
  end
  R = Icopy.base_ring
  if usefglm
     !Icopy.isGB && (Icopy = Singular.fglm(Icopy, :degrevlex));
  else
     !Icopy.isGB && (Icopy = Singular.std(Icopy, complete_reduction = true));
  end

  if Singular.gens(Icopy) == Singular.gens(oneideal(R))
    return [oneideal(R), oneideal(R)]
  end

  if Singular.vdim(Icopy) == Singular.total_degree(Icopy[1]) #Easy case, see documentation for vdim_internal
    return vdim_internal(Icopy, usefglm)
  else
    result = Vector{Singular.sideal}(undef, 0)
    return zerodecomp_internal(Icopy, result, outputReduced, usefglm)
  end
end


function zerodecomp_internal(I::Singular.sideal, result::Array{Singular.sideal}, outputReduced::Bool, usefglm::Bool)
    #Put I in general position by coordinate change with a random vector B to apply PrimaryTest, B is given back so that 
    #it can be utilized to construct the inverse map of the coordinate change
    copyI = Singular.deepcopy(I);     
    J, B = randomGeneralPosition(copyI);
    
    if usefglm
       J = Singular.fglm(J, :degrevlex);
    else
       J = Singular.std(J, complete_reduction = true);
    end
    
    R = J.base_ring;

    #g is the factorisation of the element with the smallest leading monomial
    g = Singular.factor(J[1]);
    #gi are all the (irreducible) factors of g (without multiplicities)
    gi = collect(keys(g.fac))

    #The arrays Qdash and Pdash are used to compute the elements of the primary decomposition. As general position is needed to 
    #compute and we do not want the ideal in general position for the output, P and Q denote the ideals of Pdash and Qdash obtained by applying 
    #the inverse coordinate change.
    Q = Vector{Singular.sideal}(undef, length(gi));
    Qdash = Vector{Singular.sideal}(undef, length(gi));
    P = Vector{Singular.sideal}(undef, length(gi));
    Pdash = Vector{Singular.sideal}(undef, length(gi));

    #For each factor of g another primary ideal is added to the primary decomposition (see Algorithm)
    for i in 1:length(gi)

        Qdash[i] = addGenerator(J, gi[i]^g[gi[i]]);
        Gi = Singular.Ideal(R, gi[i]);
        Q[i] = addGenerators(copyI, inverseGeneralPosition(Gi, B)^g[gi[i]]);
        Pdash[i] = primaryTest(Qdash[i], usefglm);
        if Pdash[i] != zeroideal(R)
            
            #the standard case
            P[i] = inverseGeneralPosition(Pdash[i], B);
            if outputReduced && usefglm
                result = push!(result, Singular.fglm(Q[i], :degrevlex));
                result = push!(result, Singular.fglm(P[i], :degrevlex));
            elseif outputReduced && !usefglm
                result = push!(result, Singular.std(Q[i], complete_reduction = true));
                result = push!(result, Singular.std(P[i], complete_reduction = true));
            else
                result = push!(result, Q[i]);
                result = push!(result, P[i]);
            end 
        else
            #in this case, the ideal Pdash was not in general position, i.e. we got unlucky choosing our random coordinate change.
            #repeat algorithm with only Q[i] to hopefully choose a better random coordinate change.
            result = zerodecomp_internal(Q[i], result);
        end
    end
    return result
end


################################################ PrimaryTest ############################################################################################

# This is an implementation of Algorithm 4.2.5 in "A SINGULAR introduction to Commutative Algebra" by Greuel and Pfister, 
# based somewhat on the SINGULAR implementation, after applying the preprocessing step to calculate the first generator, calls primaryTestInternal,
# which is the main function for this step

#REQUIRES:  INPUT sideal I
#           I is zero-dimensional (Krull dimension)  
#           I has a field as coefficient ring
#           basering is not a quotient ring
#           basering has a global ordering
#           basering has lex ordering
#           basering has characteristic 0
#OUTPUT:    zeroideal of basering if I is not primary or not in general position
#           radical of I otherwise
function primaryTest(I::Singular.sideal, usefglm::Bool = false)
    R = I.base_ring;
    Icopy = deepcopy(I)

    if usefglm
     !Icopy.isGB && (Icopy = Singular.fglm(Icopy, :degrevlex));
    else
     !Icopy.isGB && (Icopy = Singular.std(Icopy, complete_reduction = true));
    end

    g = Singular.factor(Icopy[1]);
    gn = collect(keys(g.fac))
    #TODO: It might be useful to implement an irreducibility check here
    for i in 1:length(gn)
      if gn[1] != gn[i]
        return Singular.Ideal(R, R(0));
      end
    end
    if length(gn) == 0
      return oneideal(R)
    end
    return primaryTestInternal(Icopy, gn[1])  
  end
  
  
  #INPUT: sideal, polynomial obtained from preprocessing
  #OUTPUT: Output for primaryTest
  function primaryTestInternal(I::Singular.sideal, p::Singular.spoly)
    R = I.base_ring;
    prm::Singular.sideal = Singular.Ideal(R, p)

    prm = Singular.std(prm; complete_reduction = true)
    n = Singular.nvars(R);
    m = 1;
    if I[1] == R(1)
      return Singular.Ideal(R, R(1))
    end
  
    while n > 1
      n-=1;  
      #Find a generating element which has a power of the n-th variable as leading term
      #Similar to the previous Singular implementation, multiply with a factor to that there appear no fractions
      m = indexOfVarN(I, n);
      e = Singular.total_degree(Singular.leading_term(I[m])); 
      t = Singular.leading_coefficient(I[m])*e*gen(R, n)+ div((Singular.tail(I[m])), gen(R, n)^(e-1));
      I[m] = (R(e)^e)*R((leading_coefficient(I[m])^(e-1)))*I[m];
      #This implies I was not in general position
      if Singular.reduce(I[m] - t^e, prm) != 0 
        return Ideal(R, R(0))
      end
     prm = addGenerator(prm, t);
  
      #reduce the primary ideal by calculating a reduced groebner basis
      prm = Singular.std(prm; complete_reduction = true)
  
    end
  
    return prm
  end
  
###################################################### helpprocs PrimaryTest ############################################################################


   # adds a generator to an ideal
  function addGenerator(I::Singular.sideal, f::Singular.spoly)
    G = Singular.gens(I);
    push!(G, f);
    J = Singular.Ideal(I.base_ring, G);
    return J
  end
  
  #adds multiple generators to an ideal
  function addGenerators(I::Singular.sideal, J::Singular.sideal)
    return I+J
  end
  
  #finds the last occurrence of a variable
  #RAUL: This is the idiot proof version, not taking into account that we have a GB wrt lex
  function indexOfVarN(I::Singular.sideal, n::Int64)
     m = 1
     while div( Singular.leading_term(I[m]),
		Singular.gen(I.base_ring, n)^(Singular.total_degree(leading_term(I[m])))) == 0
       m+=1;
     end
     return m
  end    
  
################################################# old helpprocs, should probably be deleted ###############################################################
#returns the trivial ideal
function oneideal(R::Singular.PolyRing)
    onepoly = one(R);
    res = Singular.Ideal(R, onepoly);
    return res
end

#returns the zero ideal of a PolyRing
function zeroideal(R::Singular.PolyRing)
    zeropoly = zero(R);
    res = Singular.Ideal(R, zeropoly);
    return res
end

############################################## helpprocs for coordinate changes (especially of orderings) ################################################

#changes order of basering to lex via an algebra homomorphism
function changeOrderOfBasering(I::Singular.sideal, ordering::Symbol = lex)
    R = I.base_ring;
    G = Singular.gens(I.base_ring);
    Gstrich = string.(G);
    S, G = Singular.PolynomialRing(R.base_ring, Gstrich, ordering = ordering)
    res = Singular.AlgebraHomomorphism(R, S, G);
    return res(I)
end




#requires: A has elements of the coefficient field of the base ring of the ideal of length nvars-1
function generalPosition(I::Singular.sideal, A::Vector{T}) where T <: Singular.spoly
    R = I.base_ring;
    n = Singular.nvars(R);
    G = Singular.gens(R);
    xn = Singular.gen(R, n);
    for i in 1:length(A)
        xn = xn + A[i]*gen(R, i);
    end

    G[n] = xn;

    res = Singular.AlgebraHomomorphism(R,R,G);
    return res(I)

end

#gives a random coordinate change to obtain general position (see GP, 4.2.7), calls generalPosition with 
#randomly generated elements of the basering between -5:-1, 1:5
function randomGeneralPosition(I::Singular.sideal)
    n = Singular.nvars(I.base_ring);
    A = rand(1:20, n-1)
    sign = rand(0:1, n-1)
    for i in 1:length(A)
        if sign[i] == 1
            A[i] = -A[i];
        end
    end
    B = Vector{Singular.spoly}(undef, n-1)
    R = I.base_ring;
    for i in 1:length(A)
        if A[i] == 0
            A[i] = 1
        end    
        B[i] = R(A[i]);
    end
    return [generalPosition(I,B) , B]
end

#gives the image of the ideal I under the inverse of the coordinate change obtained by general position and array A
function inverseGeneralPosition(I::Singular.sideal, A::Vector{T}) where T<: Singular.spoly
    R = I.base_ring;
    n = Singular.nvars(R);
    G = Singular.gens(R);
    xn = gen(R, n);
    for i in 1:length(A)
        xn = xn - A[i]*gen(R, i);
    end

    G[n] = xn;

    res = Singular.AlgebraHomomorphism(R,R,G);
    return res(I)

end

#############################################

##If vdim(I) == totaldegree(I[1])), I is a primary ideal and hence its decomposition is itself
function vdim_internal(I::Singular.sideal, usefglm::Bool = false)
  R = I.base_ring
  g = Singular.factor(I[1])
  gi = collect(keys(g.fac))
  primary = Vector{Singular.sideal}(undef, 0)

  #Compute ideal I without the first element for easier computation
  I_without_I1 = Singular.Ideal(R, gens(I)[2:length(gens(I))])
  for k in 1:length(gi)
       I_without_I1 = Singular.deepcopy(I_without_I1)
       prm = addGenerator(I_without_I1, gi[k]^g[gi[k]])
       p = addGenerator(I_without_I1, gi[k])
       if usefglm
          primary = push!(primary, Singular.fglm(prm, :degrevlex))
          primary = push!(primary, Singular.fglm(p, :degrevlex))
       else
          primary = push!(primary, Singular.std(prm, complete_reduction = true))
          primary = push!(primary, Singular.std(p, complete_reduction = true))
       end
  end
  return primary
end
