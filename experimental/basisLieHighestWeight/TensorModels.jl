# ?

using Oscar
#using SparseArrays

# TODO: make the first one a symmetric product, or reduce more generally

spid(n) = spdiagm(0 => [ZZ(1) for _ in 1:n])
sz(A) = size(A)[1]
tensorProduct(A, B) = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
tensorProducts(As, Bs) = (AB->tensorProduct(AB[1], AB[2])).(zip(As, Bs))
tensorPower(A, n) = (n == 1) ? A : tensorProduct(tensorPower(A, n-1), A)
tensorPowers(As, n) = (A->tensorPower(A, n)).(As)


function tensorMatricesForOperators(L, hw, ops)
    """
    Calculates the matrices g_i corresponding to the operator ops[i].
    """
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