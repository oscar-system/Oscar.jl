 


module MBOld

export basisLieHighestWeight # use parallel = true for parallel computing

using Oscar


G = Oscar.GAP.Globals
forGap = Oscar.GAP.julia_to_gap
fromGap = Oscar.GAP.gap_to_julia

#### vector space bases

TVec = SRow{ZZRingElem} 
Short = UInt8 

struct VSBasis
    A::Vector{TVec} 
    pivot::Vector{Int} 
end


nullSpace() = VSBasis([], [])


function normalize(v::TVec)
    """
    divides vector by gcd of nonzero entries, returns vector and first nonzero index
    used: addAndReduce!
    """
    if isempty(v)
        return (0, 0)
    end

    pivot = first(v)[1]

    return divexact(v, gcd(map(y->y[2], union(v)))), pivot
end


reduceCol(a, b, i::Int) = b[i]*a - a[i]*b


function addAndReduce!(sp::VSBasis, v::TVec)
    """
    for each pivot of sp.A we make entry of v zero and return the result
    0 => linear dependent
    * => linear independent, new column element of sp.A since it increases basis
    invariants: the row of a pivotelement in any column in A is 0 (except the pivotelement)
               elements of A are integers, gcd of each column is 1
    """
    A = sp.A
    pivot = sp.pivot
    v, newPivot = normalize(v) 
    if newPivot == 0
        return 0
    end
 
    for j = 1:length(A)
        i = pivot[j]
        if i != newPivot
            continue
        end
        v = reduceCol(v, A[j], i)
        v, newPivot = normalize(v)
        if newPivot == 0
            return 0
        end
    end

    pos = findfirst(pivot .> newPivot)
    if (pos === nothing)
        pos = length(pivot) + 1
    end

    insert!(A, pos, v)
    insert!(pivot, pos, newPivot)

    return v
end


#### Lie algebras

G = Oscar.GAP.Globals
forGap = Oscar.GAP.julia_to_gap
fromGap = Oscar.GAP.gap_to_julia


function lieAlgebra(t::String, n::Int)
    L = G.SimpleLieAlgebra(forGap(t), n, G.Rationals)
    return L, G.ChevalleyBasis(L)
end

# temporary workaround for issue 2128
function multiply_scalar(A::SMat{T}, d) where T
    for i in 1:nrows(A)
        scale_row!(A, i, T(d))
    end
    return A
    #return identity_matrix(SMat, QQ, size(A)[1])*A
end

gapReshape(A) = sparse_matrix(QQ, hcat(A...))

function matricesForOperators(L, hw, ops)
    """
    used to create tensorMatricesForOperators
    """
    M = G.HighestWeightModule(L, forGap(hw))
    mats = G.List(ops, o -> G.MatrixOfAction(G.Basis(M), o))
    mats = gapReshape.(fromGap(mats))
    denominators = map(y->denominator(y[2]), union(union(mats...)...))
    #d = convert(QQ, lcm(denominators))
    d = lcm(denominators)# // 1
    mats = (A->change_base_ring(ZZ, multiply_scalar(A, d))).(mats)
    return mats
end

function weightsForOperators(L, cartan, ops)
    cartan = fromGap(cartan, recursive=false)
    ops = fromGap(ops, recursive=false)
    asVec(v) = fromGap(G.ExtRepOfObj(v))
    if any(iszero.(asVec.(ops)))
        error("ops should be non-zero")
    end
    nzi(v) = findfirst(asVec(v) .!= 0)
    return [
        [asVec(h*v)[nzi(v)] / asVec(v)[nzi(v)] for h in cartan] for v in ops
    ]
end

#### tensor model

function kron(A, B)
    res = sparse_matrix(ZZ, nrows(A)*nrows(B), ncols(A)*ncols(B))
    for i in 1:nrows(B)
        for j in 1:nrows(A)
            new_row_tuples = Vector{Tuple{Int, ZZRingElem}}([(1,ZZ(0))])
            for (index_A, element_A) in union(getindex(A, j))
                for (index_B, element_B) in union(getindex(B, i))
                    push!(new_row_tuples, ((index_A-1)*ncols(B)+index_B, element_A*element_B))
                end
            end
            new_row = sparse_row(ZZ, new_row_tuples)
            setindex!(res, new_row, (j-1)*nrows(B)+i)
        end
    end
    #println("ncols(res): ", ncols(res))
    #println("nrows(res): ", nrows(res))
    
    return res
end

# temprary fix sparse in Oscar does not work
function tensorProduct(A, B)
    temp_mat = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
    res = sparse_matrix(ZZ, nrows(A)*nrows(B), ncols(A)*ncols(B))
    for i in 1:nrows(temp_mat)
        setindex!(res, getindex(temp_mat, i), i)
    end
    return res
end

# TODO: make the first one a symmetric product, or reduce more generally

#spid(n) = spdiagm(0 => [ZZ(1) for _ in 1:n])
spid(n) = identity_matrix(SMat, ZZ, n)
sz(A) = nrows(A) #size(A)[1]
#tensorProduct(A, B) = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
tensorProducts(As, Bs) = (AB->tensorProduct(AB[1], AB[2])).(zip(As, Bs))
tensorPower(A, n) = (n == 1) ? A : tensorProduct(tensorPower(A, n-1), A)
tensorPowers(As, n) = (A->tensorPower(A, n)).(As)


function tensorMatricesForOperators(L, hw, ops)
    """
    Calculates the matrices g_i corresponding to the operator ops[i].
    """
    #println("hw: ", hw)
    mats = []

    for i in 1:length(hw)
        #println("hw[i]: ", hw[i])
        if hw[i] <= 0
            continue
        end
        wi = Int.(1:length(hw) .== i) # i-th fundamental weight
        _mats = matricesForOperators(L, wi, ops)
        _mats = tensorPowers(_mats, hw[i])
        if size(mats)[1] > 0
            A = _mats[1]
            B = mats[1]
            #println("Addition cols ", ncols(kron(A, spid(sz(B))) + kron(spid(sz(A)), B)))
            #println("Addition rows ",nrows(kron(A, spid(sz(B))) + kron(spid(sz(A)), B)))
        end
        mats = mats == [] ? _mats : tensorProducts(mats, _mats)
        #println(spdiagm(0 => [ZZ(1) for _ in 1:5])agm)
        #display(mats)
    end
    return mats
end

#### monomial basis


# TODO: Demazure modules


"""
    basisLieHighestWeight(t::String, n::Int, hw::Vector{Int}; parallel::Bool = true) :: Tuple{Vector{Vector{Short}},Vector{TVec}}

Compute a monomial basis for the highest weight module with highest weight ``hw`` (in terms of the fundamental weights), for a simple Lie algebra of type ``t`` and rank ``n``.

Example
======
```jldoctest
julia> dim, monomials, vectors = PolyBases.MB.basisLieHighestWeight("A", 2, [1,0])
(3, Vector{UInt8}[[0x00, 0x00, 0x00], [0x01, 0x00, 0x00], [0x00, 0x00, 0x01]], SparseArrays.SparseVector{Int64, Int64}[  [1]  =  1,   [2]  =  1,   [3]  =  1])
```
"""

function basisLieHighestWeight(t::String, n::Int, hw::Vector{Int}; roots = []) #--- :: Tuple{Int64,Vector{Vector{Short}},Vector{TVec}}
    # TODO: flag for root ordering
    L, CH = lieAlgebra(t, n)    #--- L ist Lie Algebra vom Typ t, n (z.B: "A", 2) durch Funktion GAP.SimpleLieAlgebra(t,n)
                                #--- G ist ChevalleyBasis von L durch GAP.ChevalleyBasis(L)
                                #--- CH[1] obere strikte Dreiecksmatrizen, CH[2] untere strikte Dreiecksmatrizen, CH[3] Diagonalmatrizen

    ops = CH[1] # positive root vectors
    # .. reorder..
    wts = weightsForOperators(L, CH[3], ops)
    wts = (v->Int.(v)).(wts)
    

    # Keep only elements of ops which weights are listed in roots, in order of roots
    #if !isempty(roots)
    #    ops = fromGap(ops, recursive = false)
    #    ops_vec = []
    #    for i in 1:length(wts)
    #        if wts[i] in roots
    #            push!(ops_vec, ops[i])
    #        end
    #    end
    #    ops = forGap(ops_vec, recursive = false)
    #end 

    #--- mats speichert die Matrizen g_i für (g_1^a_1 * ... * g_k^a_k)*v0 ab.
    mats = tensorMatricesForOperators(L, hw, ops)
    #  display(wts)
    #  rnd = rand(length(wts[1]))
    #  rnd /= norm(rnd)
    #  wts = (v->round(dot(rnd,v), sigdigits=4)).(wts)
    #  display(wts)

    #display(mats)
    #display(d)
    hwv = sparse_row(ZZ, [(1,1)])
    # TODO: (okay for now)
    # here we assume the first vector in G.Basis(M) is a highest weight vector
    #display(hwv)
    #display.(1. .* (mats))

    #---
    #println("--------------------------------------")
    #println("Parameter für compute(hwv, mats, wts)")
    #println("")
    #println("hwv:")
    #display(collect(hwv))
    #println("")
    #println("mats:")
    #for i=1:length(mats)
    #    println(string("mats[", i, "]:"))
    #    #display(collect(mats[i]))
    #    display(mats[i])
    #end
    #println("")
    #println("wts:")
    #display(wts)
    #println("")

    res = compute(hwv, mats, wts)

    return length(res[1]), res...
end

nullMon(m) = zeros(Short, m)


function compute(v0, mats, wts::Vector{Vector{Int}})
    m = length(mats) #---
    monomials = [nullMon(m)]    #--- hier speichern wir die Monombasis. Ein Eintrag steht für Potenzen unserer Basiselemente
                                #--- nullMon ist leeres Monom, also nur m-Nullen bzw. die Identität
    lastPos = 0 #--- speichert für die untere while-Schleife die Anzahl der Basiselemente nach dem letzten Durchlauf

    #--- jedes Monom erhaelt id
    #--- id(mon) returnt die id des entsprechenden mon.
    id(mon) = sum((1 << (sum(mon[1:i])+i-1) for i in 1:m-1) , init = 1) # init = 1 for case m = 1
    #e = [sparse_row(ZZ, [(i,1)]) for i=1:m] #--- Einheitsvektoren z.B.: e = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] für m = 3
    e = [Short.(1:m .== i) for i in 1:m]
    #maxid(deg) = id(deg.*e[1]) #--- returnt maximale id für Monom vom Grad deg (worst case ist [1, 0, 0], da dann erste Potenz in jedem Summanden ist)
    maxid(deg) = id(deg.*e[1])
    #maxid(deg) = id(sparse_row(ZZ, [(1,deg)]))
    # TODO 

    #--- Sperrung von Monome
    blacklists = [falses(maxid(0))] #--- Monome die bereits in der Basis enthalten sind, i-ter Eintrag ist Sperrung für Monome von Grad i
    blacklists[end][id(monomials[1])] = 1 #--- Für Grad 0, d.h. konstant, ist die letzte id gesperrt da 0-Polynom
    newMons(deg) = (push!(blacklists, falses(maxid(deg)))) # Neuer Eintrag um Sperrung der Monome von Grad deg zu speichern
    checkMon(mon) = blacklists[end-1][id(mon)]  #--- returnt ob Monom gesperrt ist. Man nimmt sich den vorletzten Eintrag, weil in der Berechnung 
                                                #--- im while-Loop der letzte Eintrag immer der Eintrag für den derzeitigen Grad ist
    setMon(mon) = (blacklists[end][id(mon)] = true) #--- sperrt monom (im while loop ist deg immer maximal für mon und damit im letzten Eintrag von blacklists)


    #--- Speicher für while-loop
    vectors = [v0]  #--- von Monombasis mit (g_1^a_1 * ... * g_k^a_k)*v0 erzeugte Basisvektoren von V=R^(n+1)
    weights = [0 * wts[1]]
    space = Dict(weights[1] => nullSpace()) #--- Erzeugnis der bisherigen Basisvektoren.
    addAndReduce!(space[weights[1]], v0) #--- Fügt Span von Basisvektor v0 ein

    deg = 0

    #--- iteriere durch alle Monome in reverse lexicographical order
    #--- Bsp: [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]
    #--- while-loop iteriert durch Grad deg, also Summe der Einträge
    #--- i iteriert durch ersten nicht-null Index
    #--- di ist neue Potenz an Stelle i (d.h. 1 höher als vorher)
    #--- p durchläuft alle Monomome mit deg-1 (monomials[startPos:newPos]), um durch Addition eines e[i] alle Monome von deg erreichen
    #--- innerster Durchlauf: Erhöhe bei allen bereits verwendeten Monomen Potenz i um 1 auf di, ansonsten mache nichts
    #--- solange bis für einen Grad kein Polynom mehr hinzugefügt worden ist
    while length(monomials) > lastPos

        startPos = lastPos + 1 
        newPos = length(monomials)
        deg = deg + 1
        newMons(deg) #--- Schaffung neues Speicherplatzes in blacklists. Muss iterativ erweitert werden, weil der maximale Grad der Basis unbekannt ist

        for i in 1:m, di in deg:-1:1
            for p in startPos:newPos #--- monomials[startPos:newPos] sind alle hinzugefügten Monome von deg-1
                
                #--- Siebe Monome für rev-lex Ordnung
                if !all(monomials[p][1:i-1] .== 0) #--- Monome mit nicht-null vor Index i wurden bereits mit kleinerem i getroffen
                    continue
                end
                #--- neues monom soll an Stelle i Potenz di haben
                if monomials[p][i]+1 > di 
                    startPos = p+1 #--- wegen rev-lex Ordnung wird dieses Monom nie wieder erweitert -> für Grad deg immer erst danach beginnen
                    continue
                end
                if monomials[p][i]+1 < di
                    break #--- break, da monomials bereits rev-lex geordnet ist und alle späteren von Grad deg-1 noch kleinere Potenz bei i haben
                end

                mon = monomials[p] + e[i] #--- mögliches neues Monom

                #--- nur mon behalten, falls jeder Index != 0 durch verwendetes Monom des vorherigen Schrittes mit +e[j] hervorging.
                #--- folgt nicht automatisch induktiv, da mon-e[j] zwar überprüft worden ist, aber vielleicht neuer Vektor linear abhängig war
                if any(i != j && mon[j] > 0 && !checkMon(mon-e[j]) for j in 1:m)
                    continue
                end

                wt = weights[p] + wts[i]
                if !haskey(space, wt)
                    space[wt] = nullSpace()
                end

                vec = mul(vectors[p], transpose(mats[i])) #--- vec ist neuer potenzieller Kandidat für Basis. 
                vec = addAndReduce!(space[wt], vec) #--- vec ist altes vec ohne Anteil des bereits vorhandenen Erzeugnisses. Wenn 0, dann linear abhängig zu altem Unterraum
                if vec == 0 #--- wenn vec linear abhängig war, fügen wir vec nicht hinzu
                    continue
                end

                #display(v)
                #display(1. .* space.A)
                setMon(mon) #--- wenn vec linear unabhängig war wird vec zur Basis hinzugefügt und das zugehörige Monom gespeichert
                push!(monomials, mon)
                push!(weights, wt)
                push!(vectors, vec)
            end
        end
        lastPos = newPos
    end

    return monomials, vectors
end

end # module
