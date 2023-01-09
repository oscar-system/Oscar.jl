export depth, fitting_ideal, is_flat, is_regular_sequence, koszul_homology, non_flat_locus 

###############################################################################
# functionality supporting computations in homological algebra.
###############################################################################


##############################################################################
#
# Fitting ideals
#
##############################################################################

@doc Markdown.doc"""
     fitting_ideal(M::ModuleFP, i::Int) 

Return the `i`-th Fitting ideal of `M`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> F = free_module(R, 2);

julia> o = zero(R)
0

julia> U = matrix([x^3-y^2 o; o x^3-y^2; -x^2 y; -y x])
[x^3 - y^2           0]
[        0   x^3 - y^2]
[     -x^2           y]
[       -y           x]

julia> M = quo(F,U)[1]
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 4 generators
1 -> (x^3 - y^2)*e[1]
2 -> (x^3 - y^2)*e[2]
3 -> -x^2*e[1] + y*e[2]
4 -> -y*e[1] + x*e[2]

julia> fitting_ideal(M, -1)
ideal(0)

julia> fitting_ideal(M, 0)
ideal(x^3 - y^2)

julia> fitting_ideal(M, 1)
ideal(y, x)

julia> fitting_ideal(M, 2)
ideal(1)
```
"""
function fitting_ideal(M::ModuleFP, i::Int)
 error("not implemented for the given type of module.")
end

function fitting_ideal(F:: FreeMod, i::Int)
 R = base_ring(F)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 if i < rank(F) return ideal(R, [zero(R)]) end
 return R
end

function fitting_ideal(M::SubQuo, i::Int)
 R = base_ring(M)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 C = present_as_cokernel(M)
 U = C.quo
 SG = singular_generators(U.gens)
 return MPolyIdeal(R, Singular.LibHomolog.fitting(SG, i)) #TODO result is standard_basis
end


##############################################################################
#
# Flatness
#
##############################################################################

@doc Markdown.doc"""
     is_flat(M::ModuleFP) 

Return `true` if `M` is flat, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> F = free_module(R, 2);

julia> o = zero(R)
0

julia> U = matrix([x^3-y^2 o; o x^3-y^2; -x^2 y; -y x])
[x^3 - y^2           0]
[        0   x^3 - y^2]
[     -x^2           y]
[       -y           x]

julia> M = quo(F,U)[1]
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 4 generators
1 -> (x^3 - y^2)*e[1]
2 -> (x^3 - y^2)*e[2]
3 -> -x^2*e[1] + y*e[2]
4 -> -y*e[1] + x*e[2]

julia> is_flat(M)
false
```
"""
function is_flat(M::ModuleFP)
 error("not implemented for the given type of module.")
end

function is_flat(F:: FreeMod)
 R = base_ring(F)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 return true
end

function is_flat(M::SubQuo)
 R = base_ring(M)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 C = present_as_cokernel(M)
 U = C.quo
 SG = singular_generators(U.gens)
 res = Singular.LibHomolog.isFlat(SG)
 if res == 1 return true end
 return false
end

@doc Markdown.doc"""
     non_flat_locus(M::ModuleFP) 

Return an ideal of `base_ring(M)` which defines the non-flat-locus of `M`
in the sense that the localization of `M` at a prime ideal of `base_ring(M)`
is non-flat iff the prime ideal contains the returned ideal.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> F = free_module(R, 2);

julia> o = zero(R)
0

julia> U = matrix([x^3-y^2 o; o x^3-y^2; -x^2 y; -y x])
[x^3 - y^2           0]
[        0   x^3 - y^2]
[     -x^2           y]
[       -y           x]

julia> M = quo(F,U)[1]
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 4 generators
1 -> (x^3 - y^2)*e[1]
2 -> (x^3 - y^2)*e[2]
3 -> -x^2*e[1] + y*e[2]
4 -> -y*e[1] + x*e[2]

julia> non_flat_locus(M)
ideal(x^3 - y^2)
```
"""
function non_flat_locus(M::ModuleFP)
 error("not implemented for the given type of module.")
end

function non_flat_locus(F:: FreeMod)
 R = base_ring(F)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 return ideal(R, [one(R)])
end

function non_flat_locus(M::SubQuo)
 R = base_ring(M)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 C = present_as_cokernel(M)
 U = C.quo
 SG = singular_generators(U.gens)
 return MPolyIdeal(R, Singular.LibHomolog.flatLocus(SG))
end

##############################################################################
#
# Regular sequence test
#
##############################################################################

@doc Markdown.doc"""
     is_regular_sequence(V::Vector{T}, M::ModuleFP{T}) where T <: MPolyElem

Return `true` if the elements of `V` form, in the given order, a regular sequence on `M`.
Return `false`, otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> F = free_module(R, 1);

julia> V =  [x*z-z, x*y-y, x]
3-element Vector{fmpq_mpoly}:
 x*z - z
 x*y - y
 x

julia> is_regular_sequence(V, F)
false

julia> W = [x*z-z, x, x*y-y]
3-element Vector{fmpq_mpoly}:
 x*z - z
 x
 x*y - y

julia> is_regular_sequence(W, F)
true
```
"""
function is_regular_sequence(V::Vector{T}, M::ModuleFP{T}) where T <: MPolyElem
 error("not implemented for the given type of module.")
end

function is_regular_sequence(V::Vector{T}, F:: FreeMod{T}) where T <: MPolyElem
 R = base_ring(F)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 @assert parent(V[1]) == R
 @assert all(x->parent(x) == R, V)
 I = ideal(R, V)
 singular_assure(I)
 MX = singular_module(F)
 SG = Singular.Module(base_ring(MX), MX(zero(F)))
 res = Singular.LibHomolog.isReg(I.gens.S, SG)
 if res == 1 return true end
 return false
end

function is_regular_sequence(V::Vector{T}, M::SubQuo{T}) where T <: MPolyElem
 R = base_ring(M)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 @assert parent(V[1]) == R
 @assert all(x->parent(x) == R, V)
 I = ideal(R, V)
 singular_assure(I)
 C = present_as_cokernel(M)
 U = C.quo
 SG = singular_generators(U.gens) 
 res = Singular.LibHomolog.isReg(I.gens.S, SG)
 if res == 1 return true end
 return false
end

##############################################################################
#
# Koszul homology
#
##############################################################################

@doc Markdown.doc"""
     koszul_homology(V::Vector{T}, M::ModuleFP{T}, p::Int) where T <: MPolyElem

If $f_1, \dots, f_r$ are the entries of `V` in the given order, return the `p`-th homology 
module of the complex $K(f_1, \dots, f_r)\otimes_R M$, where $K(f_1, \dots, f_r)$ is the
*Koszul complex* defined by $f_1, \dots, f_r$.

!!! note
    See [GP08](@cite) or [DL06](@cite) for the definition of the Koszul complex.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> F = free_module(R, 1);

julia> V =  [x*y, x*z, y*z]
3-element Vector{fmpq_mpoly}:
 x*y
 x*z
 y*z

julia> koszul_homology(V, F, 0)
Submodule with 3 generators
1 -> y*z*e[1]
2 -> x*z*e[1]
3 -> x*y*e[1]
represented as subquotient with no relations.

julia> koszul_homology(V, F, 1)
Submodule with 3 generators
1 -> -z*e[1] + z*e[2]
2 -> y*e[2]
3 -> x*e[1]
represented as subquotient with no relations.

julia> koszul_homology(V, F, 2)
Submodule with 1 generator
1 -> 0
represented as subquotient with no relations.
```

```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> TC = ideal(R, [x*z-y^2, w*z-x*y, w*y-x^2]);

julia> F = free_module(R, 1);

julia> koszul_homology(gens(TC), F, 0)
Submodule with 3 generators
1 -> (-x*z + y^2)*e[1]
2 -> (-w*z + x*y)*e[1]
3 -> (-w*y + x^2)*e[1]
represented as subquotient with no relations.

julia> koszul_homology(gens(TC), F, 1)
Submodule with 3 generators
1 -> y*e[1] - z*e[2]
2 -> x*e[1] - y*e[2]
3 -> w*e[1] - x*e[2]
represented as subquotient with no relations.

julia> koszul_homology(gens(TC), F, 2)
Submodule with 1 generator
1 -> 0
represented as subquotient with no relations.
```
"""
function koszul_homology(V::Vector{T}, M::ModuleFP{T}, i::Int) where T <: MPolyElem
 error("not implemented for the given type of module.")
end

function koszul_homology(V::Vector{T}, F:: FreeMod{T}, i::Int) where T <: MPolyElem
 R = base_ring(F)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 @assert parent(V[1]) == R
 @assert all(x->parent(x) == R, V)
 I = ideal(R, V)
 singular_assure(I)
 MX = singular_module(F)
 SG = Singular.Module(base_ring(MX), MX(zero(F)))
 SU = Singular.LibHomolog.KoszulHomology(I.gens.S, SG, i)
 FFF = free_module(R, rank(SU))
 MG = ModuleGens(FFF, SU)
 UO = SubModuleOfFreeModule(FFF, MG)
 return SubQuo(UO)
end

function koszul_homology(V::Vector{T}, M::SubQuo{T}, i::Int) where T <: MPolyElem
 R = base_ring(M)
 if !(coefficient_ring(R) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 @assert parent(V[1]) == R
 @assert all(x->parent(x) == R, V)
 I = ideal(R, V)
 singular_assure(I)
 C = present_as_cokernel(M)
 U = C.quo
 SG = singular_generators(U.gens)
 SU = Singular.LibHomolog.KoszulHomology(I.gens.S, SG, i)
 FFF = free_module(R, rank(SU))
 MG = ModuleGens(FFF, SU)
 UO = SubModuleOfFreeModule(FFF, MG)
 return SubQuo(UO)
end

##############################################################################
#
# depth
#
##############################################################################

@doc Markdown.doc"""
     depth(I::MPolyIdeal{T}, M::ModuleFP{T}) where T <: MPolyElem

Return the depth of `I` on `M`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> TC = ideal(R, [x*z-y^2, w*z-x*y, w*y-x^2]);

julia> dim(TC)
2

julia> F = free_module(R, 1);

julia> U = collect(gen(TC, i)*F[1] for i in 1:ngens(TC));

julia> M, _ = quo(F, U);

julia> I = ideal(R, gens(R))
ideal(w, x, y, z)

julia> depth(I, M)
2
```

```jldoctest
julia> S, x, y = PolynomialRing(QQ, "x" => 1:3, "y" => 1:5);

julia> W = [y[1]-x[1]^2,  y[2]-x[2]^2,   y[3]-x[3]^2, y[4]-x[2]*(x[1]-x[3]),  y[5]-(x[1]-x[2])*x[3]];

julia> J = eliminate(ideal(S, W), x);

julia> R, y = PolynomialRing(QQ, "y" => 1:5);

julia> W = append!(repeat([zero(R)], 3), gens(R))
8-element Vector{fmpq_mpoly}:
 0
 0
 0
 y[1]
 y[2]
 y[3]
 y[4]
 y[5]

julia> P = hom(S, R, W);

julia> VP4 = P(J);

julia> dim(VP4)
3

julia> F = free_module(R, 1);

julia> U = collect(gen(VP4, i)*F[1] for i in 1:ngens(VP4));

julia> M, _ = quo(F, U);

julia> I = ideal(R, gens(R))
ideal(y[1], y[2], y[3], y[4], y[5])

julia> depth(I, M)
1
```
"""
function depth(I::MPolyIdeal{T}, M::ModuleFP{T}) where T <: MPolyElem
 error("not implemented for the given type of module.")
end

function depth(I::MPolyIdeal{T}, F:: FreeMod{T}) where T <: MPolyElem
 if !(coefficient_ring(base_ring(I)) isa AbstractAlgebra.Field)
    throw(ArgumentError("The coefficient ring must be a field."))
 end
 @assert base_ring(I) == base_ring(F)
 if is_zero(F) || is_equal((I*F)[1], F) return -1 end
 MX = singular_module(F)
 SG = Singular.Module(base_ring(MX), MX(zero(F)))
 singular_assure(I)
 return Singular.LibHomolog.depth(SG, I.gens.S)
end

function depth(I::MPolyIdeal{T}, M::SubQuo{T}) where T <: MPolyElem
 if !(coefficient_ring(base_ring(I)) isa AbstractAlgebra.Field)
   throw(ArgumentError("The coefficient ring must be a field."))
 end
 @assert base_ring(I) == base_ring(M)
 if is_zero(M) || is_equal((I*M)[1], M) return -1 end
 singular_assure(I)
 C = present_as_cokernel(M)
 U = C.quo
 SG = singular_generators(U.gens)
 return Singular.LibHomolog.depth(SG, I.gens.S)
end


