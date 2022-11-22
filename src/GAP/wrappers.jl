# Use the @wrap macro to set up optimized versions of a bunch of frequently
# used GAP functions. So instead of writing e.g. GAP.Globals.DenominatorRat(x)
# you'd write GAPWrap.DenominatorRat(x). The former always performs a variable
# lookup in the GAP kernel and is not type stable. The latter accesses the
# underlying GAP function object directly and also has information about the
# return type.
#
# Note that the macro GAP.@wrap has a similar purpose as @gapwrap and @gappatribute
# have, but works on a much lower level on purpose. We may actually phase out
# use of @gapwrap and @gappatribute in the future.
module GAPWrap

using GAP

GAP.@wrap CF(x::Any, y::Any)::GapObj
GAP.@wrap CF(x::Any)::GapObj
GAP.@wrap CHAR_FFE_DEFAULT(x::Any)::GapInt
GAP.@wrap Characteristic(x::Any)::GapInt
GAP.@wrap Coefficients(x::Any, y::Any)::GapObj
GAP.@wrap Conductor(x::Any)::GapInt
GAP.@wrap CycList(x::GapObj)::GapInt
GAP.@wrap DegreeFFE(x::Any)::Int
GAP.@wrap DenominatorCyc(x::Any)::GapInt
GAP.@wrap DenominatorRat(x::Any)::GapInt
GAP.@wrap E(x::Any)::GapInt
GAP.@wrap ExtRepOfObj(x::GapObj)::GapObj
GAP.@wrap GF(x::Any, y::Any)::GapObj
GAP.@wrap GF(x::Any)::GapObj
GAP.@wrap Image(x::Any, y::Any)::GapObj
GAP.@wrap Image(x::Any)::GapObj
GAP.@wrap IndependentGeneratorExponents(x::Any, y::Any)::GapObj
GAP.@wrap INT_FFE_DEFAULT(x::Any)::GapInt
GAP.@wrap IntFFE(x::Any)::GapInt
GAP.@wrap IsAbelian(x::Any)::Bool
GAP.@wrap IsAlgebraicExtension(x::Any)::Bool
GAP.@wrap IsAlmostSimpleGroup(x::Any)::Bool
GAP.@wrap IsAlternatingForm(x::Any)::Bool
GAP.@wrap IsAlternatingGroup(x::Any)::Bool
GAP.@wrap IsBiCoset(x::Any)::Bool
GAP.@wrap IsBijective(x::Any)::Bool
GAP.@wrap IsCharacteristicSubgroup(x::Any, y::Any)::Bool
GAP.@wrap IsCheapConwayPolynomial(x::Any, y::Any)::Bool
GAP.@wrap IsClassFunction(x::Any)::Bool
GAP.@wrap IsConjugate(x::Any, y::Any, z::Any)::Bool
GAP.@wrap IsCyc(x::Any)::Bool
GAP.@wrap IsCyclic(x::Any)::Bool
GAP.@wrap IsCyclotomic(x::Any)::Bool
GAP.@wrap IsCyclotomicCollColl(x::Any)::Bool
GAP.@wrap IsDihedralGroup(x::Any)::Bool
GAP.@wrap IsDoneIterator(x::Any)::Bool
GAP.@wrap IsEmpty(x::Any)::Bool
GAP.@wrap IsFFE(x::Any)::Bool
GAP.@wrap IsFinite(x::Any)::Bool
GAP.@wrap IsGroupOfAutomorphisms(x::Any)::Bool
GAP.@wrap IsHermitianForm(x::Any)::Bool
GAP.@wrap IsInjective(x::Any)::Bool
GAP.@wrap IsInnerAutomorphism(x::Any)::Bool
GAP.@wrap IsInt(x::Any)::Bool
GAP.@wrap IsList(x::Any)::Bool
GAP.@wrap IsMatrix(x::GapObj)::Bool
GAP.@wrap IsMatrixGroup(x::GapObj)::Bool
GAP.@wrap IsMatrixObj(x::GapObj)::Bool
GAP.@wrap IsMatrixOrMatrixObj(x::Any)::Bool
GAP.@wrap IsNaturalAlternatingGroup(x::Any)::Bool
GAP.@wrap IsNaturalSymmetricGroup(x::Any)::Bool
GAP.@wrap IsNilpotentGroup(x::Any)::Bool
GAP.@wrap IsNormal(x::Any, y::Any)::Bool
GAP.@wrap IsOne(x::Any)::Bool
GAP.@wrap IsPcGroup(x::Any)::Bool
GAP.@wrap IsPerfectGroup(x::Any)::Bool
GAP.@wrap IsPermGroup(x::Any)::Bool
GAP.@wrap IsPGroup(x::Any)::Bool
GAP.@wrap IsPrimitive(x::Any, y::Any)::Bool
GAP.@wrap IsPrimitive(x::Any)::Bool
GAP.@wrap IsQuaternionGroup(x::Any)::Bool
GAP.@wrap IsRegular(x::Any, y::Any)::Bool
GAP.@wrap IsRegular(x::Any)::Bool
GAP.@wrap IsSemiRegular(x::Any, y::Any)::Bool
GAP.@wrap IsSemiRegular(x::Any)::Bool
GAP.@wrap IsSimpleGroup(x::Any)::Bool
GAP.@wrap IsSingularForm(x::Any)::Bool
GAP.@wrap IsSolvableGroup(x::Any)::Bool
GAP.@wrap IsSubgroupFpGroup(x::Any)::Bool
GAP.@wrap IsSubset(x::Any, y::Any)::Bool
GAP.@wrap IsSupersolvableGroup(x::Any)::Bool
GAP.@wrap IsSurjective(x::Any)::Bool
GAP.@wrap IsSymmetricForm(x::Any)::Bool
GAP.@wrap IsSymmetricGroup(x::Any)::Bool
GAP.@wrap IsTransitive(x::Any, y::Any)::Bool
GAP.@wrap IsTransitive(x::Any)::Bool
GAP.@wrap IsZero(x::Any)::Bool
GAP.@wrap LargestMovedPoint(x::Any)::Int
GAP.@wrap NextIterator(x::GapObj)::Any
GAP.@wrap NrCols(x::GapObj)::Int
GAP.@wrap NrConjugacyClasses(x::Any)::GapInt
GAP.@wrap NrRows(x::GapObj)::Int
GAP.@wrap NumberColumns(x::GapObj)::Int
GAP.@wrap NumberRows(x::GapObj)::Int
GAP.@wrap NumeratorRat(x::Any)::GapInt
GAP.@wrap One(x::Any)::Any
GAP.@wrap Order(x::Any)::GapInt
GAP.@wrap Size(x::Any)::GapInt
GAP.@wrap UnderlyingElement(x::GapObj)::GapObj
GAP.@wrap Zero(x::Any)::Any

# for Int arguments we can sometimes provide better alternatives
Conductor(x::Int) = 1
DenominatorCyc(x::Int) = 1
DenominatorRat(x::Int) = 1
IsFFE(x::Int) = false
IsOne(x::Int) = isone(x)
IsZero(x::Int) = iszero(x)
NumeratorRat(x::Int) = x
One(x::Int) = 1
Zero(x::Int) = 0

end
