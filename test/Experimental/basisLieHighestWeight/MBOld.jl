#--- Bens "Lernkommentare" mit #--- oder """...""" in Methode gekenzeichnet 


module MBOld

export basisLieHighestWeight # use parallel = true for parallel computing

using Oscar
using SparseArrays

#### vector space bases

ZZ = Int
TVec = SparseVector{ZZ, Int} #--- Werte sind ZZ, Indizies sind Int (TVec ist Datentyp der Basisvektoren in basisLieHighestWeight)
Short = UInt8 # for exponents of monomials; max. 255

struct VSBasis
    A::Vector{TVec} #--- Vektor von Basisvektoren
    pivot::Vector{Int} #--- Vektor von Pivotelementen, d.h. pivot[i] ist erster nichtnull Index von A[i]
end


nullSpace() = VSBasis([], []) #--- leerer Vektorraum


function normalize(v::TVec)
    """
    divides vector by gcd of nonzero entries, returns vector and first nonzero index
    used: addAndReduce!
    """
    dropzeros!(v) #--- deletes stored zeros in SparsedArray (in-place, hingegen würde dropzeros(v) Kopie zurückgeben)
    if isempty(v.nzind) #--- v.nzind sind abgespeicherte Werte
                        #--- v.nzval sind abgespeicherte Indizies (der Größe nach geordnet, startet bei 1)
        return (0, 0)
    end

    pivot = v.nzind[1]  #--- erster != 0 Eintrag im Vektor, d.h. der oberste

    return v .÷ gcd(v.nzval), pivot #--- durch ggT teilen
end


reduceCol(a, b, i::Int) = a*b[i] - b*a[i] #--- a = reduceCol(a, b, i) erzeugt Nulleintrag in a[i] durch elementare Zeilenumformungen
                                          #--- Werte bleiben Integers


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
    if newPivot == 0 #--- v ist Nullvektor
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


gapReshape(A) = sparse(hcat(A...))


function matricesForOperators(L, hw, ops)
    M = G.HighestWeightModule(L, forGap(hw))
    mats = G.List(ops, o -> G.MatrixOfAction(G.Basis(M), o))
    mats = gapReshape.(fromGap(mats))
    d = lcm(denominator.(union(mats...)))
    mats = (A->ZZ.(A*d)).(mats)
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

# TODO: make the first one a symmetric product, or reduce more generally

spid(n) = spdiagm(0 => [ZZ(1) for _ in 1:n])
sz(A) = size(A)[1]
tensorProduct(A, B) = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
tensorProducts(As, Bs) = (AB->tensorProduct(AB[1], AB[2])).(zip(As, Bs))
tensorPower(A, n) = (n == 1) ? A : tensorProduct(tensorPower(A, n-1), A)
tensorPowers(As, n) = (A->tensorPower(A, n)).(As)


function tensorMatricesForOperators(L, hw, ops)
    mats = []

    for i in 1:length(hw)
        if hw[i] <= 0
            continue
        end
        wi = Int.(1:length(hw) .== i) # i-th fundamental weight
        _mats = matricesForOperators(L, wi, ops)
        _mats = tensorPowers(_mats, hw[i])
        mats = mats == [] ? _mats : tensorProducts(mats, _mats)
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

function basisLieHighestWeight(t::String, n::Int, hw::Vector{Int}; roots = [], parallel::Bool = false) #--- :: Tuple{Int64,Vector{Vector{Short}},Vector{TVec}}
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

    d = sz(mats[1])
    #display(mats)
    #display(d)
    hwv = spzeros(ZZ, d); hwv[1] = 1
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

    res = compute(hwv, mats, wts, parallel = parallel)

    return length(res[1]), res...
end

nullMon(m) = zeros(Short, m)


######### compute_0 ###################
# serial computing
function compute_0(v0, mats, wts::Vector{Vector{Int}})
    m = length(mats) #---
    monomials = [nullMon(m)]    #--- hier speichern wir die Monombasis. Ein Eintrag steht für Potenzen unserer Basiselemente
                                #--- nullMon ist leeres Monom, also nur m-Nullen bzw. die Identität
    lastPos = 0 #--- speichert für die untere while-Schleife die Anzahl der Basiselemente nach dem letzten Durchlauf

    #--- jedes Monom erhaelt id
    #--- id(mon) returnt die id des entsprechenden mon.
    id(mon) = sum((1 << (sum(mon[1:i])+i-1) for i in 1:m-1) , init = 1) # init = 1 for case m = 1
    e = [Short.(1:m .== i) for i in 1:m] #--- Einheitsvektoren z.B.: e = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] für m = 3
    maxid(deg) = id(deg.*e[1]) #--- returnt maximale id für Monom vom Grad deg (worst case ist [1, 0, 0], da dann erste Potenz in jedem Summanden ist)
    
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

                vec = mats[i] * vectors[p] #--- vec ist neuer potenzieller Kandidat für Basis. 
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


############ compute_p ########################
# parallel computing
function compute_p(v0, mats, wts::Vector{Vector{Int}})
    m = length(mats)
    monomials = [nullMon(m)]
    lastPos = 0

    id(mon) = sum((1 << (sum(mon[1:i])+i-1) for i in 1:m-1), init = 1) # init = 1 for case m = 1
    e = [Short.(1:m .== i) for i in 1:m]
    maxid(deg) = id(deg.*e[1])

    blacklists = [falses(maxid(0))]
    blacklists[end][id(monomials[1])] = 1
    newMons(deg) = (push!(blacklists, falses(maxid(deg))))
    checkMon(mon) = blacklists[end-1][id(mon)]
    setMon(mon) = (blacklists[end][id(mon)] = true)

    vectors = [v0]
    weights = [0 * wts[1]]
    idweight = Dict(weights[1] => 1)
    nweights = 1
    space = Dict(1 => nullSpace())
    addAndReduce!(space[1], v0)

    deg = 0
    while length(monomials) > lastPos
        startPos = lastPos+1
        newPos = length(monomials)
        deg = deg + 1
        newMons(deg)

        num = 0
        jobs = Dict{Int, Vector{Tuple{Int, Vector{Short}, Vector{Int64}, SparseMatrixCSC{Int, Int}, SparseVector{Int, Int}, VSBasis}}}()

        for i in 1:m, di in deg:-1:1
            for p in startPos:newPos
                # lex iterator
                if !all(monomials[p][1:i-1] .== 0) ||
                        monomials[p][i] + 1 > di
                        startPos = p + 1
                    continue
                elseif monomials[p][i] + 1 < di
                    break
                end

                mon = monomials[p] + e[i]

                if any(i != j && mon[j] > 0 && !checkMon(mon - e[j]) for j in 1:m)
                    continue
                end

                num += 1
                wt = weights[p] + wts[i]

                if !haskey(idweight, wt)
                    nweights += 1
                    idweight[wt] = nweights
                end
                idwt = idweight[wt]

                if !haskey(space, idwt)
                    space[idwt] = nullSpace()
                end

                jo = (num, mon, wt, mats[i], vectors[p], space[idwt])
                if !haskey(jobs, idwt)
                    jobs[idwt] = [jo]
                else
                    push!(jobs[idwt], jo)
                end
                #display(typeof(jo))
            end
        end
        #display(typeof(jobs))

        #display(Threads.nthreads())

        res = Vector{Any}(nothing, num)
        idwts = [k for (k, l) in jobs if !isempty(l)]
        Threads.@threads for idwt in idwts
            for (i, mon, wt, A, v, sp) in jobs[idwt]
                # this is the heavy lifting
                # (in particular, the matrix product A * v)
                w = addAndReduce!(sp, A * v)
                if w != 0
                    res[i] = (mon, wt, w)
                end
            end
        end

        #display(res)

        for r in res
            if r == nothing
                continue
            end
            (mon, wt, vec) = r
            setMon(mon)
            push!(monomials, mon)
            push!(vectors, vec)
            push!(weights, wt)
        end

        lastPos = newPos
    end

    #display(idweight)
    #display(vectors)
    #display([(i, sp.A) for (i, sp) in space if length(sp.pivot) > 0])
    return monomials, vectors
end


##################
# chooses parallel / serial computing
# parallel true to use parallel computing through function compute_p 
compute(v0, mats, wts::Vector{Vector{Int}}; parallel::Bool = false) = parallel ? compute_p(v0, mats, wts) : compute_0(v0, mats, wts)  

end # module
