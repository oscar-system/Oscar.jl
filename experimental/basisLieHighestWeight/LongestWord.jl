#using Gapjm

G = Oscar.GAP.Globals
forGap = Oscar.GAP.julia_to_gap
fromGap = Oscar.GAP.gap_to_julia

#function longest_weyl_word(t, n)
#    """
#    generates a reduced expression of the longest weylword of type (t,n) by choosing uniformly a random reflection that is not leftdescending
#    the resulting longest words are not uniformly distributed
#    """
#    W = Gapjm.coxgroup(Symbol(t),n) # Weyl-group
#    S = Gapjm.gens(W) # generators of W (simple reflections)
#    m = length(S)
#    p = W() # id
#    word = [] # generated word
#
#    # extend p with reflection that are not leftdescending until not possible (i.e. reached longest word)
#    while true
#        not_desc = [i for i=1:m if !(i in Gapjm.leftdescents(W,p))] # set of i s.t. length(S[i]*p) > length(p)
#        if length(not_desc) >= 1
#            i = rand(not_desc)
#            push!(word, i)
#            p = S[i]*p
#        else
#            break
#        end
#    end
#    return word
#end

#function is_longest_weyl_word(t,n,word)
#    """
#    returns if word is a reduced expression of the longest weyl word of type (t,n) 
#    is_longest_weyl_word(t,n,longest_weyl_word(t, n)) is always true
#    """
#    W = Gapjm.coxgroup(Symbol(t),n) # Weyl-group
#    p = W(word ...) # group element of word
#    return p == longest(W) # is word longest word?
#end

function sub_simple_refl(word, L, n)
    """
    substitute simple reflections (i,i+1), saved in dec by i, with E_{i,i+1}  
    """
    R = G.RootSystem(L)
    CG = fromGap(G.CanonicalGenerators(R)[1], recursive = false)
    ops = forGap([CG[i] for i in word], recursive = false)
    return ops
end
