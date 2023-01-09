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

GAP.@wrap AsList(x::GapObj)::GapObj
GAP.@wrap AsSet(x::GapObj)::GapObj
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
GAP.@wrap Elements(x::GapObj)::GapObj
GAP.@wrap ExtRepOfObj(x::GapObj)::GapObj
GAP.@wrap FreeAbelianGroup(x::Int)::GapObj
GAP.@wrap FreeGeneratorsOfFpGroup(x::GapObj)::GapObj
GAP.@wrap FreeGroupOfFpGroup(x::GapObj)::GapObj
GAP.@wrap GeneratorsOfGroup(x::GapObj)::GapObj
GAP.@wrap GF(x::Any, y::Any)::GapObj
GAP.@wrap GF(x::Any)::GapObj
GAP.@wrap GroupHomomorphismByFunction(x1, x2, x3)::GapObj
GAP.@wrap GroupHomomorphismByFunction(x1, x2, x3, x4)::GapObj
GAP.@wrap GroupHomomorphismByFunction(x1, x2, x3, x4, x5)::GapObj
GAP.@wrap Image(x::Any, y::Any)::GapObj
GAP.@wrap Image(x::Any)::GapObj
GAP.@wrap IndependentGeneratorExponents(x::Any, y::Any)::GapObj
GAP.@wrap Inverse(x::GapObj)::GapObj
GAP.@wrap INT_FFE_DEFAULT(x::Any)::GapInt
GAP.@wrap IntFFE(x::Any)::GapInt
GAP.@wrap IsAbelian(x::Any)::Bool
GAP.@wrap IsAlgebraicElementCollCollColl(x::Any)::Bool
GAP.@wrap IsAlgebraicExtension(x::Any)::Bool
GAP.@wrap IsAlmostSimpleGroup(x::Any)::Bool
GAP.@wrap IsAlternatingForm(x::Any)::Bool
GAP.@wrap IsAlternatingGroup(x::Any)::Bool
GAP.@wrap IsAssocWord(x::Any)::Bool
GAP.@wrap IsBiCoset(x::Any)::Bool
GAP.@wrap IsBijective(x::Any)::Bool
GAP.@wrap IsBool(x::Any)::Bool
GAP.@wrap IsChar(x::Any)::Bool
GAP.@wrap IsCharacteristicSubgroup(x::Any, y::Any)::Bool
GAP.@wrap IsCheapConwayPolynomial(x::Any, y::Any)::Bool
GAP.@wrap IsCheapConwayPolynomial(x::Any)::Bool
GAP.@wrap IsClassFunction(x::Any)::Bool
GAP.@wrap IsConjugate(x::Any, y::Any, z::Any)::Bool
GAP.@wrap IsCyc(x::Any)::Bool
GAP.@wrap IsCyclic(x::Any)::Bool
GAP.@wrap IsCyclotomic(x::Any)::Bool
GAP.@wrap IsCyclotomicCollColl(x::Any)::Bool
GAP.@wrap IsCyclotomicCollCollColl(x::Any)::Bool
GAP.@wrap IsCyclotomicCollection(x::Any)::Bool
GAP.@wrap IsCyclotomicField(x::Any)::Bool
GAP.@wrap IsDihedralGroup(x::Any)::Bool
GAP.@wrap IsDoneIterator(x::Any)::Bool
GAP.@wrap IsDuplicateTable(x::Any)::Bool
GAP.@wrap IsElementaryAbelian(x::Any)::Bool
GAP.@wrap IsEmpty(x::Any)::Bool
GAP.@wrap IsFFE(x::Any)::Bool
GAP.@wrap IsFFECollCollColl(x::Any)::Bool
GAP.@wrap IsField(x::Any)::Bool
GAP.@wrap IsFinite(x::Any)::Bool
GAP.@wrap IsFinitelyGeneratedGroup(x::Any)::Bool
GAP.@wrap IsFreeGroup(x::Any)::Bool
GAP.@wrap IsFpGroup(x::Any)::Bool
GAP.@wrap IsGroupOfAutomorphisms(x::Any)::Bool
GAP.@wrap IsHandledByNiceMonomorphism(x::Any)::Bool
GAP.@wrap IsHermitianForm(x::Any)::Bool
GAP.@wrap IsInjective(x::Any)::Bool
GAP.@wrap IsInnerAutomorphism(x::Any)::Bool
GAP.@wrap IsInt(x::Any)::Bool
GAP.@wrap IsIntegers(x::Any)::Bool
GAP.@wrap IsIrreducibleCharacter(x::Any)::Bool
GAP.@wrap IsLetterAssocWordRep(x::Any)::Bool
GAP.@wrap IsLetterWordsFamily(x::Any)::Bool
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
GAP.@wrap IsPolynomial(x::Any)::Bool
GAP.@wrap IsPrimeField(x::Any)::Bool
GAP.@wrap IsPrimitive(x::Any, y::Any)::Bool
GAP.@wrap IsPrimitive(x::Any)::Bool
GAP.@wrap IsQuasisimpleGroup(x::Any)::Bool
GAP.@wrap IsQuaternionGroup(x::Any)::Bool
GAP.@wrap IsRationals(x::Any)::Bool
GAP.@wrap IsRegular(x::Any, y::Any)::Bool
GAP.@wrap IsRegular(x::Any)::Bool
GAP.@wrap IsSemiRegular(x::Any, y::Any)::Bool
GAP.@wrap IsSemiRegular(x::Any)::Bool
GAP.@wrap IsSet(x::Any)::Bool
GAP.@wrap IsSimpleGroup(x::Any)::Bool
GAP.@wrap IsSingularForm(x::Any)::Bool
GAP.@wrap IsSolvableGroup(x::Any)::Bool
GAP.@wrap IsSporadicSimpleGroup(x::Any)::Bool
GAP.@wrap IsSubgroupFpGroup(x::Any)::Bool
GAP.@wrap IsSubset(x::Any, y::Any)::Bool
GAP.@wrap IsSupersolvableGroup(x::Any)::Bool
GAP.@wrap IsSurjective(x::Any)::Bool
GAP.@wrap IsSyllableAssocWordRep(x::Any)::Bool
GAP.@wrap IsSyllableWordsFamily(x::Any)::Bool
GAP.@wrap IsSymmetricForm(x::Any)::Bool
GAP.@wrap IsSymmetricGroup(x::Any)::Bool
GAP.@wrap IsTransitive(x::Any, y::Any)::Bool
GAP.@wrap IsTransitive(x::Any)::Bool
GAP.@wrap IsTrivial(x::Any)::Bool
GAP.@wrap IsUnivariatePolynomialRing(x::Any)::Bool
GAP.@wrap IsWholeFamily(x::Any)::Bool
GAP.@wrap IsZero(x::Any)::Bool
GAP.@wrap IsZmodnZObj(x::Any)::Bool
GAP.@wrap IsZmodnZObjNonprimeCollection(x::Any)::Bool
GAP.@wrap Iterator(x::Any)::GapObj
GAP.@wrap LargestMovedPoint(x::Any)::Int
GAP.@wrap NextIterator(x::GapObj)::Any
GAP.@wrap NrCols(x::GapObj)::Int
GAP.@wrap NrConjugacyClasses(x::Any)::GapInt
GAP.@wrap NrRows(x::GapObj)::Int
GAP.@wrap NumberColumns(x::GapObj)::Int
GAP.@wrap NumberRows(x::GapObj)::Int
GAP.@wrap NumeratorRat(x::Any)::GapInt
GAP.@wrap ObjByExtRep(x::GapObj, y::GapObj)::GapObj
GAP.@wrap One(x::Any)::GAP.Obj
GAP.@wrap Order(x::Any)::GapInt
GAP.@wrap PrimePGroup(x::GapObj)::GapInt
GAP.@wrap Range(x::GapObj)::GapObj
GAP.@wrap RelatorsOfFpGroup(x::GapObj)::GapObj
GAP.@wrap Size(x::Any)::GapInt
GAP.@wrap Source(x::GapObj)::GapObj
GAP.@wrap StringViewObj(x::Any)::GapObj
GAP.@wrap UnderlyingElement(x::GapObj)::GapObj
GAP.@wrap Z(x::Any)::GAP.Obj
GAP.@wrap Zero(x::Any)::GAP.Obj

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
