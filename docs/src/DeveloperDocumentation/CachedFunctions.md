# Functions Which Cache Things

Many functions take in a keyword argument `cached::Bool`, which then cache or do not cache something which I am not fully qualified to comment on. (FIXME). This can sometimes lead to strange behaviour. So we document all the functions which do this, and at a future date, maybe do something with this list. 

## Oscar.jl

28 Functions

> - `FreeMod(R::Ring, n::Int, name::VarName = :e; cached::Bool = false)`
> - `FreeMod(R::Ring, names::Vector{String}; cached::Bool=false)`
> - `FreeMod(R::Ring, names::Vector{Symbol}; cached::Bool=false)`
> - `FreeMod_dec(R::CRing_dec, n::Int, name::VarName = :e; cached::Bool = false)`
> - `FreeMod_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::VarName = :e; cached::Bool = false)`
> - `kaehler_differentials(R::Union{MPolyRing, MPolyLocRing}; cached::Bool=true)`
> - `kaehler_differentials(R::MPolyQuoRing; cached::Bool=true)`
> - `kaehler_differentials(R::MPolyQuoLocRing; cached::Bool=true)`
> - `kaehler_differentials(R::Ring, p::Int; cached::Bool=true)`
> - `de_rham_complex(R::Ring; cached::Bool=true)`
> - `exterior_power(M::SubquoModule, p::Int; cached::Bool=true)`
> - `exterior_power(F::FreeMod, p::Int; cached::Bool=true)`
> - `koszul_complex(v::FreeModElem; cached::Bool=true)`
> - `koszul_complex(v::FreeModElem, M::ModuleFP; cached::Bool=true)`
> - `koszul_homology(v::FreeModElem, i::Int; cached::Bool=true)`
> - `koszul_homology(v::FreeModElem, M::ModuleFP, i::Int; cached::Bool=true)`
> - `polynomial_ring(R::TropicalSemiring, s::Symbol; cached::Bool = true)`
> - `polynomial_ring(R::AbstractAlgebra.Ring, v1::Pair{<:VarName, <:Any}, v...; cached::Bool = false, ordering::Symbol = :lex)`
> - `Hecke.number_field(::QQField, a::QQAbElem; cached::Bool = false)`
> - `Hecke.number_field(::QQField, a::AbstractVector{<: QQAbElem}; cached::Bool = false)`
> - `SLPolyRing(r::Ring, s::Vector{Symbol}; cached::Bool = false)`
> - `Hecke.number_field(::QQField, chi::GAPGroupClassFunction; cached::Bool = false)`
> - `slpoly_ring(R::AbstractAlgebra.Ring, n::Int; cached::Bool = false)`
> - `slpoly_ring(R::AbstractAlgebra.Ring, p::Pair{Symbol, <:AbstractVector{Int}}...; cached::Bool = false)`
> - `extension_field(f::ZZPolyRingElem, n::String = "_a"; cached::Bool = true, check::Bool = true)`
> - `extension_field(f::QQPolyRingElem, n::String = "_a"; cached::Bool = true, check::Bool = true)`
> - `extension_field(f::Generic.Poly{<:Generic.RationalFunctionFieldElem{T}}, n::String = "_a";  cached::Bool = true, check::Bool = true)`
> - `extension_field(f::Generic.Poly{nf_elem}, n::String = "_a";  cached::Bool = true, check::Bool = true)`

## Hecke.jl

80 Functions

> - `GroupAlgebra(K::Ring, G::FinGenAbGroup, cached::Bool = true)`
> - `Order(A::S, B::Vector{T}; check::Bool = true, isbasis::Bool = false, cached::Bool = true)`
> - `Order(A::S, M::FakeFmpqMat; check::Bool = true, cached::Bool = true)`
> - `MaximalOrder(O::AlgAssAbsOrd{S, T}; cached::Bool = true)`
> - `function_field(f::PolyRingElem{<:Generic.RationalFunctionFieldElem}, s::VarName = :_a; check::Bool = true, cached::Bool = false)`
> - `extension_field(f::PolyRingElem{<:Generic.RationalFunctionFieldElem}, s::VarName = :_a; check::Bool = true, cached::Bool = false)`
> - `prime_decomposition(O::GenOrd, p::RingElem, degree_limit::Int = degree(O)`
> - `FacElemMon{S}(R::S, cached::Bool = false)`
> - `FacElemMon{AbsSimpleNumField}(R::AbsSimpleNumField, cached::Bool = true)`
> - `AbsNumFieldOrderSet{T}(a::T, cached::Bool = false)`
> - `AbsNumFieldOrder{S, T}(K::S, x::FakeFmpqMat, xinv::FakeFmpqMat, B::Vector{T}, cached::Bool = false)`
> - `AbsNumFieldOrder{S, T}(K::S, x::FakeFmpqMat, cached::Bool = false)`
> - `AbsNumFieldOrder{S, T}(b::Vector{T}, cached::Bool = false)`
> - `AbsNonSimpleNumField(ff::Vector{QQPolyRingElem}, f::Vector{QQMPolyRingElem}, S::Vector{Symbol}, cached::Bool = false)`
> - `KInftyRing{T}(K::Generic.RationalFunctionField{T}, cached::Bool)`
> - `eisenstein_extension(f::Generic.Poly{S}, s::String = "a"; check::Bool = true, cached::Bool = true)`
> - `unramified_extension(f::Generic.Poly{S}, s::String = "a"; check::Bool = true, cached::Bool = true)`
> - `local_field(f::Generic.Poly{S},::Type{T}; check::Bool = true, cached::Bool = true)`
> - `local_field(f::Generic.Poly{S}, s::String, ::Type{EisensteinLocalField}; check::Bool = true, cached::Bool = true)`
> - `local_field(f::Generic.Poly{S}, s::String, ::Type{UnramifiedLocalField}; check::Bool = true, cached::Bool = true)`
> - `local_field(f::Generic.Poly{S}, s::String, ::Type{T} = GenericLocalField; check::Bool = true, cached::Bool = true)`
> - `local_field(f::QQPolyRingElem, p::Int, precision::Int, s::String, ::Type{T} = GenericLocalField; check::Bool = true, cached::Bool = true)`
> - `image(mF::NfToFqNmodMor_easy, a::FacElem{AbsSimpleNumFieldElem, AbsSimpleNumField}, D::Vector, cached::Bool, quo::Int = 0)`
> - `image(mF::NfToGFMor_easy, a::FacElem{AbsSimpleNumFieldElem, AbsSimpleNumField}, D::Vector, cached::Bool, quo::Int = 0)`
> - `image(mF::NfToGFMor_easy, a::AbsSimpleNumFieldElem, D::Vector, cached::Bool, n_quo::Int = 0)`
> - `OrdLoc{T}(OK::AbsNumFieldOrder{AbsSimpleNumField,T}, prime::AbsNumFieldOrderIdeal{AbsSimpleNumField,T}, cached::Bool = true, comp::Bool = false)`
> - `number_field(S::EuclideanRingResidueRing{QQPolyRingElem}; cached::Bool = true, check::Bool = true)`
> - `number_field(f::ZZPolyRingElem, s::Symbol; cached::Bool = true, check::Bool = true)`
> - `number_field(f::ZZPolyRingElem, s::AbstractString; cached::Bool = true, check::Bool = true)`
> - `number_field(f::ZZPolyRingElem; cached::Bool = true, check::Bool = true)`
> - `radical_extension(n::Int, gen::Integer; cached::Bool = true, check::Bool = true)`
> - `radical_extension(n::Int, gen::ZZRingElem; cached::Bool = true, check::Bool = true)`
> - `wildanger_field(n::Int, B::ZZRingElem, s::String = "_\$"; check::Bool = true, cached::Bool = true)`
> - `wildanger_field(n::Int, B::Integer, s::String = "_\$"; cached::Bool = true, check::Bool = true)`
> - `quadratic_field(d::IntegerUnion; cached::Bool = true, check::Bool = true)`
> - `quadratic_field(d::ZZRingElem; cached::Bool = true, check::Bool = true)`
> - `quadratic_field(d::Integer; cached::Bool = true, check::Bool = true)`
> - `simple_extension(K::AbsNonSimpleNumField; cached::Bool = true, check = true, simplified::Bool = false)`
> - `number_field(K1::AbsSimpleNumField, K2::AbsSimpleNumField; cached::Bool = false, check::Bool = false)`
> - `number_field(fields::Vector{AbsSimpleNumField}; cached::Bool = true, check::Bool = true)`
> - `number_field(f::Vector{QQPolyRingElem}, s::String="_\$"; cached::Bool = false, check::Bool = true)`
> - `number_field(f::Vector{QQPolyRingElem}, s::Vector{String}; cached::Bool = false, check::Bool = true)`
> - `number_field(f::Vector{QQPolyRingElem}, S::Vector{Symbol}; cached::Bool = false, check::Bool = true)`
> - `number_field(f::Vector{ZZPolyRingElem}, s::String="_\$"; cached::Bool = false, check::Bool = true)`
> - `number_field(f::Vector{ZZPolyRingElem}, s::Vector{String}; cached::Bool = false, check::Bool = true)`
> - `number_field(f::Vector{ZZPolyRingElem}, S::Vector{Symbol}; cached::Bool = false, check::Bool = true)`
> - `cyclotomic_field(::Type{NonSimpleNumField}, n::Int, s::String="z"; cached::Bool = false)`
> - `simplify(K::AbsSimpleNumField; canonical::Bool = false, cached::Bool = true, save_LLL_basis::Bool = true)`
> - `number_field(f::PolyRingElem{<: NumFieldElem}; cached::Bool = false, check::Bool = true)`
> - `number_field(::Type{AbsSimpleNumField}, L::RelSimpleNumField{AbsSimpleNumFieldElem}; check::Bool = true, cached::Bool = true)`
> - `cyclotomic_field_as_cm_extension(n::Int; cached::Bool = true)`
> - `number_field(f::Vector{Generic.Poly{T}}, S::Vector{Symbol}; cached::Bool = false, check::Bool = true)`
> - `number_field(f::Vector{Generic.Poly{T}}, s::String="_\$"; cached::Bool = false, check::Bool = true)`
> - `simplify(K::RelSimpleNumField; cached::Bool = true, prec::Int = 100)`
> - `RelNonSimpleNumField(abs_pol::Array{Generic.Poly{T}}, f::Vector{Nemo.Generic.MPoly{T}}, S::Vector{Symbol}; cached::Bool = false)`
> - `absolute_simple_field(K::AbsNonSimpleNumField; cached::Bool = true, simplify::Bool = false)`
> - `absolute_simple_field(K::NumField; cached::Bool = false, simplify::Bool = false)`
> - `absolute_simple_field(K::RelSimpleNumField{AbsSimpleNumFieldElem}; cached::Bool = false, simplify::Bool = false)`
> - `collapse_top_layer(K::RelSimpleNumField{T}; cached::Bool = false, do_embedding::Bool = true)`
> - `prime_decomposition(O::AbsNumFieldOrder{<:NumField{QQFieldElem}, <:Any}, p::IntegerUnion, degree_limit::Int = degree(O)`
> - `prime_decomposition(O::AbsSimpleNumFieldOrder, p::IntegerUnion, degree_limit::Int = degree(O)`
> - `norm_change_const(O::AbsSimpleNumFieldOrder; cached::Bool = true)`
> - `extend(O::AbsNumFieldOrder, elts::Vector{T}; check::Bool = true, cached::Bool = true)`
> - `extend(O::RelNumFieldOrder, elts::Vector{T}; check::Bool = true, cached::Bool = true)`
> - `EquationOrder(K::NumField{QQFieldElem}, cached::Bool = true)`
> - `EquationOrder(f::ZZPolyRingElem; cached::Bool = true, check::Bool = true)`
> - `EquationOrder(f::QQPolyRingElem; cached::Bool = true, check::Bool = true)`
> - `+(a::AbsNumFieldOrder, b::AbsNumFieldOrder; cached::Bool = false)`
> - `hermitian_space(E::NumField, n::Int; cached::Bool = true)`
> - `hermitian_space(E::NumField, gram::MatElem; cached::Bool = true)`
> - `rescale(V::HermSpace, r, cached::Bool=true)`
> - `quadratic_space(K::Field, n::Int; cached::Bool = true)`
> - `quadratic_space(K::Field, G::MatElem; check::Bool = true, cached::Bool = true)`
> - `rescale(q::QuadSpace, r; cached::Bool = true)`
> - `norm_group(f::Nemo.PolyRingElem, mR::T, is_abelian::Bool = true; of_closure::Bool = false, cached::Bool = true, check::Bool = false)`
> - `norm_group(l_pols::Vector{T}, mR::U, is_abelian::Bool = true; of_closure::Bool = false, cached::Bool = true, check::Bool = false)`
> - `cyclotomic_extension(k::AbsSimpleNumField, n::Int; cached::Bool = true, compute_maximal_order::Bool = true, compute_LLL_basis::Bool = true, simplified::Bool = true)`
> - `cyclotomic_extension(k::ClassField, n::Int; cached::Bool = true)`
> - `cyclotomic_extension(::Type{ClassField}, zk::AbsSimpleNumFieldOrder, n::Int; cached::Bool = true)`
> - `cyclotomic_extension(::Type{ClassField}, k::AbsSimpleNumField, n::Int; cached::Bool = true)`


## Nemo

95 Functions

> - `AbsSimpleNumField(pol::QQPolyRingElem, s::Symbol, cached::Bool = false, check::Bool = true)`
> - `number_field(f::QQPolyRingElem, s::VarName = "_a"; cached::Bool = true, check::Bool = true)`
> - `residue_field(R::QQPolyRing, f::QQPolyRingElem; cached::Bool = true)`
> - `ArbField(p::Int = 256; cached::Bool = true)`
> - `AcbField(p::Int = 256; cached::Bool = true)`
> - `RealPolyRing(R::RealField, S::Symbol, cached::Bool = true)`
> - `ArbPolyRing(R::ArbField, S::Symbol, cached::Bool = true)`
> - `ComplexPolyRing(R::ComplexField, S::Symbol, cached::Bool = true)`
> - `AcbPolyRing(R::AcbField, S::Symbol, cached::Bool = true)`
> - `similar(f::PolyRingElem, R::ComplexField, var::VarName=var(parent(f)`
> - `similar(f::PolyRingElem, R::RealField, var::VarName=var(parent(f)`
> - `similar(f::PolyRingElem, R::AcbField, var::VarName=var(parent(f)`
> - `similar(f::PolyRingElem, R::ArbField, var::VarName=var(parent(f)`
> - `similar(f::PolyRingElem, R::FpField, s::Symbol=var(parent(f)`
> - `similar(f::PolyRingElem, R::fpField, s::Symbol=var(parent(f)`
> - `residue_ring(R::ZZRing, n::ZZRingElem; cached::Bool=true)`
> - `residue_ring(R::ZZRing, n::Integer; cached::Bool=true)`
> - `finite_field(q::IntegerUnion, s::VarName = :o; cached::Bool = true, check::Bool = true)`
> - `finite_field(f::FqPolyRingElem, s::VarName = :o; cached::Bool = true, check::Bool = true, absolute::Bool = false)`
> - `GF(a::IntegerUnion, s::VarName = :o; cached::Bool = true, check::Bool = true)`
> - `GF(p::IntegerUnion, d::Int, s::VarName = :o; cached::Bool = true, check::Bool = true)`
> - `GF(f::FqPolyRingElem, s::VarName = :o; cached::Bool = true, check::Bool = true, absolute::Bool = false)`
> - `residue_field(R::ZZRing, p::IntegerUnion; cached::Bool = true)`
> - `FqField(f::FqPolyRingElem, s::Symbol, cached::Bool = false, absolute::Bool = false)`
> - `residue_ring(R::ZZRing, n::Int; cached::Bool=true)`
> - `residue_ring(R::ZZRing, n::UInt; cached::Bool=true)`
> - `QadicField(p::Integer, d::Int, prec::Int, var::String = "a"; cached::Bool = true)`
> - `ZZPolyRing(R::ZZRing, s::Symbol, cached::Bool = true)`
> - `QQPolyRing(R::QQField, s::Symbol, cached::Bool = true)`
> - `zzModRing(n::UInt, cached::Bool=true)`
> - `ZZModRing(n::ZZRingElem, cached::Bool=true)`
> - `zzModPolyRing(R::zzModRing, s::Symbol, cached::Bool = true)`
> - `fpPolyRing(R::fpField, s::Symbol, cached::Bool = true)`
> - `ZZModPolyRing(R::ZZModRing, s::Symbol, cached::Bool = true)`
> - `FpPolyRing(R::FpField, s::Symbol, cached::Bool = true)`
> - `ZZMPolyRing(s::Vector{Symbol}, S::Symbol, cached::Bool = true)`
> - `QQMPolyRing(s::Vector{Symbol}, S::Symbol, cached::Bool = true)`
> - `zzModMPolyRing(R::zzModRing, s::Vector{Symbol}, S::Symbol, cached::Bool = true)`
> - `fpMPolyRing(R::fpField, s::Vector{Symbol}, S::Symbol, cached::Bool = true)`
> - `FpMPolyRing(R::FpField, s::Vector{Symbol}, S::Symbol, cached::Bool = true)`
> - `fqPolyRepField(c::ZZRingElem, deg::Int, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `fqPolyRepField(f::zzModPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `fqPolyRepField(f::fpPolyRingElem, s::Symbol, cached::Bool = true; check::Bool=true)`
> - `FqField(char::ZZRingElem, deg::Int, s::Symbol, cached::Bool = true)`
> - `FqField(f::ZZModPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqField(f::FpPolyRingElem, s::Symbol, cached::Bool = true;  check::Bool = true)`
> - `FqField(f::zzModPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqField(f::fpPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqPolyRepField(char::ZZRingElem, deg::Int, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqPolyRepField(f::ZZModPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqPolyRepField(f::FpPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqPolyRepField(char::ZZRingElem, deg::Int, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqPolyRepField(f::ZZModPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `FqPolyRepField(f::FpPolyRingElem, s::Symbol, cached::Bool = true; check::Bool = true)`
> - `fqPolyRepMPolyRing(R::fqPolyRepField, s::Vector{Symbol}, S::Symbol, cached::Bool = true)`
> - `PadicField(p::ZZRingElem, prec::Int; cached::Bool = true, check::Bool = true)`
> - `QadicField(p::ZZRingElem, d::Int, prec::Int, var::String = "a"; cached::Bool = true, check::Bool = true)`
> - `ZZRelPowerSeriesRing(prec::Int, s::Symbol, cached::Bool = true)`
> - `ZZAbsPowerSeriesRing(prec::Int, s::Symbol, cached::Bool = true)`
> - `FlintPuiseuxSeriesRing{T}(R::Ring, cached::Bool = true)`
> - `FlintPuiseuxSeriesField{T}(R::Field, cached::Bool = true)`
> - `ZZLaurentSeriesRing(prec::Int, s::Symbol, cached::Bool = true)`
> - `QQRelPowerSeriesRing(prec::Int, s::Symbol, cached::Bool = true)`
> - `QQAbsPowerSeriesRing(prec::Int, s::Symbol, cached::Bool = true)`
> - `FqPolyRing(R::FqField, s::Symbol, cached::Bool = true)`
> - `FqPolyRepPolyRing(R::FqPolyRepField, s::Symbol, cached::Bool = true)`
> - `fqPolyRepPolyRing(R::fqPolyRepField, s::Symbol, cached::Bool = true)`
> - `FqMPolyRing(R::FqField, s::Vector{Symbol}, internal_ordering::Symbol = :lex, cached::Bool = true)`
> - `similar(f::PolyRingElem, R::QQField, s::Symbol=var(parent(f)`
> - `matrix_space(R::ZZRing, r::Int, c::Int; cached::Bool = true)`
> - `matrix_space(R::ZZModRing, r::Int, c::Int; cached::Bool = true)`
> - `similar(f::PolyRingElem, R::ZZModRing, s::Symbol=var(parent(f)`
> - `similar(f::PolyRingElem, R::ZZRing, s::Symbol=var(parent(f)`
> - `matrix_space(R::FqField, r::Int, c::Int; cached::Bool = true)`
> - `similar(f::PolyRingElem, R::FqField, s::Symbol=var(parent(f)`
> - `matrix_space(R::FqPolyRepField, r::Int, c::Int; cached::Bool = true)`
> - `matrix_space(R::fqPolyRepField, r::Int, c::Int; cached::Bool = true)`
> - `similar(f::PolyRingElem, R::fqPolyRepField, s::Symbol=var(parent(f)`
> - `similar(f::PolyRingElem, R::FqPolyRepField, s::Symbol=var(parent(f)`
> - `matrix_space(R::FpField, r::Int, c::Int; cached::Bool = true)`
> - `matrix_space(R::fpField, r::Int, c::Int; cached::Bool = true)`
> - `matrix_space(R::zzModRing, r::Int, c::Int; cached::Bool = true)`
> - `similar(f::PolyRingElem, R::zzModRing, s::Symbol=var(parent(f)`
> - `prime_field(F::fqPolyRepField; cached::Bool=true)`
> - `prime_field(F::FqPolyRepField; cached::Bool=true)`
> - `prime_field(F::T; cached::Bool=true)`
> - `GF(n::Int; cached::Bool=true, check::Bool=true)`
> - `GF(n::UInt; cached::Bool=true, check::Bool=true)`
> - `GF(n::ZZRingElem; cached::Bool=true, check::Bool=true)`
> - `finite_field(F::FqPolyRepField, deg::Int, s::VarName = :o; cached::Bool = true, check::Bool=true)`
> - `finite_field(char::Int, deg::Int, s::VarName = :o; cached::Bool = true, check::Bool = true)`
> - `finite_field(pol::Zmodn_poly, s::VarName = :o; cached::Bool = true, check::Bool=true)`
> - `finite_field(F::fqPolyRepField, deg::Int, s::VarName = :o; cached::Bool = true, check::Bool = true)`
> - `finite_field(p::Integer; cached::Bool = true, check::Bool = true)`
> - `finite_field(p::ZZRingElem; cached::Bool = true, check::Bool = true)`

## Abstract Algebra

62 Functions

> - `rational_function_field(k::Field, s::VarName; cached::Bool=true)`
> - `residue_ring(R::Ring, a::RingElement; cached::Bool = true)`
> - `residue_ring(R::PolyRing, a::RingElement; cached::Bool = true)`
> - `quo(R::Ring, a::RingElement; cached::Bool = true)`
> - `residue_field(R::Ring, a::RingElement; cached::Bool = true)`
> - `quo(::Type{Field}, R::Ring, a::RingElement; cached::Bool = true)`
> - `SparsePolynomialRing(R::Ring, s::VarName; cached::Bool = true)`
> - `total_ring_of_fractions(R::Ring; cached::Bool=true)`
> - `FactoredFractionField(R::AbstractAlgebra.Ring; cached::Bool=true)`
> - `function_field(p::Poly{RationalFunctionFieldElem{T, U}}, s::VarName; cached::Bool=true)`
> - `function_field(p::Generic.Poly{Generic.RationalFunctionFieldElem{T, U}}, s::VarName; cached::Bool=true)`
> - `map_coefficients(g::T, p::LaurentMPolyWrap; cached::Bool = true,`
> - `laurent_polynomial_ring(R::AbstractAlgebra.Ring, s::Vector{Symbol}; cached::Bool = true)`
> - `change_base_ring(R::Ring, p::LaurentPolyWrap; cached::Bool = true,`
> - `map_coefficients(g::T, p::LaurentPolyWrap; cached::Bool = true,`
> - `laurent_polynomial_ring(R::AbstractAlgebra.Ring, s::Symbol; cached::Bool = true)`
> - `laurent_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true)`
> - `laurent_series(R::Field, arr::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true)`
> - `laurent_series_ring(R::AbstractAlgebra.Ring, prec::Int, s::VarName; cached::Bool=true)`
> - `laurent_series_ring(R::AbstractAlgebra.Field, prec::Int, s::VarName; cached::Bool=true)`
> - `laurent_series_field(R::AbstractAlgebra.Field, prec::Int, s::VarName; cached::Bool=true)`
> - `LocalizedEuclideanRing{T}(prime::T, primes::Vector{T}, cached::Bool = true, comp::Bool = false)`
> - `LocalizedEuclideanRing{T}(prime::T, cached::Bool = true, comp::Bool = false)`
> - `localization(R::AbstractAlgebra.Ring, prime::T; cached::Bool=true, comp::Bool = false)`
> - `localization(R::AbstractAlgebra.Ring, primes::Vector{T}; cached::Bool=true)`
> - `power_series_ring(R::AbstractAlgebra.Ring, prec::Int, s::VarName; cached::Bool=true, model=:capped_relative)`
> - `SparsePolynomialRing(R::AbstractAlgebra.Ring, s::Symbol; cached::Bool = true)`
> - `total_ring_of_fractions(R::Ring; cached::Bool=true)`
> - `fraction_field(R::AbstractAlgebra.Ring; cached::Bool=true)`
> - `PolyRing{T}(R::Ring, s::Symbol, cached::Bool = true)`
> - `MPolyRing{T}(R::Ring, s::Vector{Symbol}, internal_ordering::Symbol=:lex, cached::Bool=true)`
> - `UniversalPolyRing{T, U}(R::Ring, ord::Symbol, cached::Bool=true)`
> - `EuclideanRingResidueRing{T}(modulus::T, cached::Bool = true)`
> - `EuclideanRingResidueField{T}(modulus::T, cached::Bool = true)`
> - `RelPowerSeriesRing{T}(R::Ring, prec::Int, s::Symbol, cached::Bool = true)`
> - `AbsPowerSeriesRing{T}(R::Ring, prec::Int, s::Symbol, cached::Bool = true)`
> - `FracField{T}(R::Ring, cached::Bool = true)`
> - `MatSpace{T}(R::NCRing, r::Int, c::Int, cached::Bool = true)`
> - `FreeAssAlgebra{T}(R::Ring, s::Vector{Symbol}, cached::Bool = true)`
> - `IdealSet{T}(R::Ring, cached::Bool = true)`
> - `matrix_ring(R::AbstractAlgebra.NCRing, n::Int; cached::Bool = true)`
> - `matrix_space(R::AbstractAlgebra.NCRing, r::Int, c::Int; cached::Bool = true)`
> - `rational_function_field(k::Field, s::Symbol; cached::Bool=true)`
> - `rational_function_field(k::Field, s::Vector{Symbol}; cached::Bool=true)`
> - `change_base_ring(R::Ring, p::UnivPoly{T, U}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p)`
> - `change_coefficient_ring(R::Ring, p::UnivPoly{T, U}; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(R, parent(p)`
> - `map_coefficients(f::T, p::UnivPoly; cached::Bool=true, parent::UniversalPolyRing = _change_univ_poly_ring(parent(f(zero(base_ring(p)`
> - `universal_polynomial_ring(R::Ring; internal_ordering=:lex, cached::Bool=true)`
> - `GF(p::T; cached::Bool = true, check::Bool=true)`
> - `abs_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true)`
> - `fraction_field(R::Ring; cached::Bool=true)`
> - `FactoredFractionField(R::Ring; cached::Bool=true)`
> - `free_module(R::NCRing, rank::Int; cached::Bool = true)`
> - `vector_space(R::Field, dim::Int; cached::Bool = true)`
> - `change_base_ring(R::Ring, p::MPolyRingElem{T}; cached::Bool=true, parent::MPolyRing = _change_mpoly_ring(R, parent(p)`
> - `change_coefficient_ring(R::Ring, p::MPolyRingElem{T}; cached::Bool=true, parent::MPolyRing = _change_mpoly_ring(R, parent(p)`
> - `change_base_ring(R::NCRing, p::NCPolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p)`
> - `change_coefficient_ring(R::NCRing, p::NCPolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p)`
> - `polynomial(R::Ring, arr::Vector{T}, var::VarName=:x; cached::Bool=true)`
> - `change_base_ring(R::Ring, p::PolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p)`
> - `change_coefficient_ring(R::Ring, p::PolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p)`
> - `rel_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true)`

## Singular

10 Functions

> - `N_ZnRing(n::BigInt, cached::Bool = true)`
> - `N_ZpField(n::Int, cached::Bool = true)`
> - `N_GField(p::Int, n::Int, S::Symbol, cached::Bool = true)`
> - `N_FField(F::Field, S::Vector{Symbol}, cached::Bool = true)`
> - `N_FField(F::Field, S::AbstractVector{<:VarName}, cached::Bool = true)`
> - `N_AlgExtField(F::N_FField, minpoly::n_transExt, cached::Bool = true)`
> - `FunctionField(F::Singular.Field, S::AbstractVector{<:VarName}; cached::Bool=true)`
> - `FunctionField(F::Singular.Field, n::Int; cached::Bool=true)`
> - `CoefficientRing(R::T, cached::Bool = true)`
> - `GAlgebra(R::PolyRing{T}, C, D; cached::Bool = true)`

## GAP

No caching here!

## Polymake

No caching here!


## Generation script

This output can be generated by running the following script inside the root folder of each project:

```bash
#!/bin/bash
grep -d recurse "^export" src | sed "s/^.*export /\^/" | sed "s/, /\n/g" | sed "s/\(\^[0-9a-zA-Z_\!]*\).*/\1/"  > exports.txt
grep -n -d recurse -w "^\s*function" src | grep "cached::Bool" | sed "s/^.*function //" | sed "s/).*$/)/" | grep -f exports.txt -w
```