# Use the @wrap macro to set up optimized versions of a bunch of frequently
# used GAP functions. So instead of writing e.g. GAP.Globals.DenominatorRat(x)
# you'd write GAPWrap.DenominatorRat(x). The former always performs a variable
# lookup in the GAP kernel and is not type stable. The latter accesses the
# underlying GAP function object directly and also has information about the
# return type.
#
# Note that the macro GAP.@wrap has a similar purpose as @gapattribute has,
# but works on a much lower level on purpose. We may actually phase out
# use of @gapattribute in the future.
module GAPWrap

using GAP

# the following list is intended to be sorted according to `LC_COLLATE=C sort -f`, i.e. ignoring case
GAP.@wrap AbelianGroup(x::GapObj, y::GapObj)::GapObj
GAP.@wrap AbelianPcpGroup(x::GAP.Obj, y::GapObj)::GapObj
GAP.@wrap AlgebraicExtension(x::GapObj, y::GapObj)::GapObj
GAP.@wrap AlgExtElm(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap AntiSymmetricParts(x::GapObj, y::GapObj, z::GapInt)::GapObj
GAP.@wrap AsList(x::GapObj)::GapObj
GAP.@wrap AsSet(x::GapObj)::GapObj
GAP.@wrap AssocWordByLetterRep(x::GapObj, y::GapObj)::GapObj
GAP.@wrap AssocWordByLetterRep(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap AtlasIrrationality(x::GapObj)::GAP.Obj
GAP.@wrap AutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap Basis(x::GapObj)::GapObj
GAP.@wrap Basis(x::GapObj, y::GapObj)::GapObj
GAP.@wrap BasisNC(x::GapObj, y::GapObj)::GapObj
GAP.@wrap BrauerCharacterValue(x::GapObj)::GAP.Obj
GAP.@wrap CanonicalBasis(x::GapObj)::GapObj
GAP.@wrap CentralCharacter(x::GapObj)::GapObj
GAP.@wrap CentreOfCharacter(x::GapObj, y::GapObj)::GapObj
GAP.@wrap CF(x::Any)::GapObj
GAP.@wrap CF(x::Any, y::Any)::GapObj
GAP.@wrap Characteristic(x::Any)::GapInt
GAP.@wrap CharacterParameters(x::GapObj)::GapObj
GAP.@wrap CharacterTable(x::GapObj)::GapObj
GAP.@wrap CharacterTable(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap CharacterTableFactorGroup(x::GapObj, y::GapObj)::GapObj
GAP.@wrap CharacterTableWreathSymmetric(x::GapObj, y::GapInt)::GapObj
GAP.@wrap CHAR_FFE_DEFAULT(x::Any)::GapInt
GAP.@wrap ClassFunction(x::GapObj, y::GapObj)::GapObj
GAP.@wrap ClassMultiplicationCoefficient(x::GapObj, y::Int, z::Int, t::Int)::GAP.Obj
GAP.@wrap ClassNames(x::GapObj)::GapObj
GAP.@wrap ClassParameters(x::GapObj)::GapObj
GAP.@wrap ClassPositionsOfCentre(x::GapObj)::GapObj
GAP.@wrap ClassPositionsOfDerivedSubgroup(x::GapObj)::GapObj
GAP.@wrap ClassPositionsOfNormalSubgroups(x::GapObj)::GapObj
GAP.@wrap ClassPositionsOfPCore(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap ClassPositionsOfSolvableResiduum(x::GapObj)::GapObj
GAP.@wrap Coefficients(x::Any, y::Any)::GapObj
GAP.@wrap CoefficientsFamily(x::GapObj)::GapObj
GAP.@wrap CoefficientsOfUnivariatePolynomial(x::GapObj)::GapObj
GAP.@wrap CoeffsCyc(x::GAP.Obj, y::Int)::GapObj
GAP.@wrap CollectionsFamily(x::GapObj)::GapObj
GAP.@wrap Collector(x::GapObj)::GapObj
GAP.@wrap ComputedClassFusions(x::GapObj)::GapObj
GAP.@wrap ComputedPowerMaps(x::GapObj)::GapObj
GAP.@wrap Conductor(x::Any)::GapInt
GAP.@wrap ConjugacyClass(x::GapObj, y::GapObj)::GapObj
GAP.@wrap ConjugacyClasses(x::GapObj)::GapObj
GAP.@wrap ConjugacyClassesMaximalSubgroups(x::GapObj)::GapObj
GAP.@wrap ConjugacyClassesSubgroups(x::GapObj)::GapObj
GAP.@wrap ConjugacyClassSubgroups(x::GapObj, y::GapObj)::GapObj
GAP.@wrap ConjugateSubgroup(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Core(x::GapObj, y::GapObj)::GapObj
GAP.@wrap CycleFromList(x::GapObj)::GapObj
GAP.@wrap CycleStructurePerm(x::GapObj)::GapObj
GAP.@wrap CYCLE_LENGTH_PERM_INT(x::GapObj, y::Int)::Int
GAP.@wrap CycList(x::GapObj)::GapInt
GAP.@wrap CyclotomicPol(x::Int)::GapObj
GAP.@wrap DecomposedFixedPointVector(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DecomposeTensorProduct(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap Decomposition(x::GapObj, y::GapObj, z::GAP.Obj)::GapObj
GAP.@wrap DecompositionMatrix(x::GapObj)::GapObj
GAP.@wrap DefiningPolynomial(x::GapObj)::GapObj
GAP.@wrap DegreeFFE(x::Any)::Int
GAP.@wrap DegreeOfLaurentPolynomial(x::GapObj)::GapInt
GAP.@wrap DegreeOverPrimeField(x::GapObj)::Int
GAP.@wrap Depth(x::GapObj)::Int
GAP.@wrap DepthOfPcElement(x::GapObj, y::GapObj)::Int
GAP.@wrap DenominatorCyc(x::Any)::GapInt
GAP.@wrap DenominatorRat(x::Any)::GapInt
GAP.@wrap DescriptionOfRootOfUnity(x::Any)::GapObj
GAP.@wrap DeterminantOfCharacter(x::GapObj)::GapObj
GAP.@wrap Dimension(x::GapObj)::Int
GAP.@wrap DimensionOfHighestWeightModule(x::GapObj, y::GapObj)::GapInt
GAP.@wrap DominantCharacter(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DoubleCoset(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap DoubleCosetRepsAndSizes(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap E(x::Any)::GapInt
GAP.@wrap EigenvaluesChar(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap ElementOfFpGroup(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Elements(x::GapObj)::GapObj
GAP.@wrap ElementsFamily(x::GapObj)::GapObj
GAP.@wrap ELMS_LIST(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Embedding(x::GapObj, y::Int)::GapObj
GAP.@wrap EpimorphismSchurCover(x::GapObj)::GapObj
GAP.@wrap Exponents(x::GapObj)::GapObj
GAP.@wrap ExponentsOfPcElement(x::GapObj, y::GapObj)::GapObj
GAP.@wrap ExtraspecialGroup(x::GapObj, y::GAP.Obj, z::GapObj)::GapObj
GAP.@wrap ExtRepOfObj(x::GapObj)::GapObj
GAP.@wrap ExtRepPolynomialRatFun(x::GapObj)::GapObj
GAP.@wrap FactorCosetAction(x::GapObj, y::GapObj)::GapObj
GAP.@wrap FamilyObj(x::GAP.Obj)::GapObj
GAP.@wrap FamilyPcgs(x::GAP.Obj)::GapObj
GAP.@wrap fhmethsel(x::GapObj)::GAP.Obj
GAP.@wrap Field(x::Any)::GapObj
GAP.@wrap Flat(x::GapObj)::GapObj
GAP.@wrap FreeAbelianGroup(x::Int)::GapObj
GAP.@wrap FreeGeneratorsOfFpGroup(x::GapObj)::GapObj
GAP.@wrap FreeGroupOfFpGroup(x::GapObj)::GapObj
GAP.@wrap FusionCharTableTom(x::GapObj, y::GapObj)::GapObj
GAP.@wrap FusionConjugacyClasses(x::GapObj, y::GapObj)::GapObj
GAP.@wrap GaloisCyc(x::GAP.Obj, GapInt)::GAP.Obj
GAP.@wrap GeneratorsOfAlgebra(x::GapObj)::GapObj
GAP.@wrap GeneratorsOfField(x::GapObj)::GapObj
GAP.@wrap GeneratorsOfGroup(x::GapObj)::GapObj
GAP.@wrap GenExpList(x::GapObj)::GapObj
GAP.@wrap GetFusionMap(x::GapObj, y::GapObj)::GapObj
GAP.@wrap GF(x::Any)::GapObj
GAP.@wrap GF(x::Any, y::Any)::GapObj
GAP.@wrap Group(x::GapObj)::GapObj
GAP.@wrap Group(x::GapObj, y::GapObj)::GapObj
GAP.@wrap GroupHomomorphismByFunction(x1, x2, x3)::GapObj
GAP.@wrap GroupHomomorphismByFunction(x1, x2, x3, x4)::GapObj
GAP.@wrap GroupHomomorphismByFunction(x1, x2, x3, x4, x5)::GapObj
GAP.@wrap GroupOfPcgs(x::GapObj)::GapObj
GAP.@wrap Grp(x::GapObj)::GapObj
GAP.@wrap HasCharacterParameters(x::GapObj)::Bool
GAP.@wrap HasClassParameters(x::GapObj)::Bool
GAP.@wrap HasConjugacyClassesSubgroups(x::GapObj)::Bool
GAP.@wrap Hasfhmethsel(x::GapObj)::Bool
GAP.@wrap HasGrp(x::GapObj)::Bool
GAP.@wrap HasImageRecogNode(x::GapObj)::Bool
GAP.@wrap HasIsRecogInfoForAlmostSimpleGroup(x::GapObj)::Bool
GAP.@wrap HasIsRecogInfoForSimpleGroup(x::GapObj)::Bool
GAP.@wrap HasKernelRecogNode(x::GapObj)::Bool
GAP.@wrap HasMaxes(x::GapObj)::Bool
GAP.@wrap HasMaximalAbelianQuotient(x::Any)::Bool
GAP.@wrap HasSize(x::Any)::Bool
GAP.@wrap HirschLength(x::GapObj)::Int
GAP.@wrap Identifier(x::GapObj)::GapObj
GAP.@wrap Identity(x::GapObj)::GapObj
GAP.@wrap Image(x::Any)::GapObj
GAP.@wrap Image(x::Any, y::Any)::GapObj
GAP.@wrap ImageRecogNode(x::GapObj)::GapObj
GAP.@wrap ImagesRepresentative(x::GapObj, y::Any)::GAP.Obj
GAP.@wrap ImagesSource(x::GapObj)::GapObj
GAP.@wrap ImmutableMatrix(x::GapObj, y::GapObj, z::Bool)::GapObj
GAP.@wrap IndependentGeneratorExponents(x::Any, y::Any)::GapObj
GAP.@wrap IndependentGeneratorsOfAbelianGroup(x::GapObj)::GapObj
GAP.@wrap Indeterminate(x::GapObj)::GapObj
GAP.@wrap Indeterminate(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap IndeterminateNumberOfUnivariateRationalFunction(x::GapObj)::Int
GAP.@wrap IndeterminatesOfPolynomialRing(x::GapObj)::GapObj
GAP.@wrap Indicator(x::GapObj, y::GapInt)::GapObj
GAP.@wrap Indicator(x::GapObj, y::GapObj, z::GapInt)::GapObj
GAP.@wrap InducedClassFunction(x::GapObj, y::GapObj)::GapObj
GAP.@wrap InducedClassFunctionsByFusionMap(x::GapObj, y::GapObj, z::GapObj, t::GapObj)::GapObj
GAP.@wrap InducedCyclic(x::GapObj)::GapObj
GAP.@wrap InducedCyclic(x::GapObj, y::GapObj)::GapObj
GAP.@wrap InitFusion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Intersection(x::GapObj)::GapObj
GAP.@wrap IntFFE(x::Any)::GapInt
GAP.@wrap INT_FFE_DEFAULT(x::Any)::GapInt
GAP.@wrap Inverse(x::GapObj)::GapObj
GAP.@wrap Irr(x::GapObj)::GapObj
GAP.@wrap IsAbelian(x::Any)::Bool
GAP.@wrap IsAlgebraicElementCollCollColl(x::Any)::Bool
GAP.@wrap IsAlgebraicExtension(x::Any)::Bool
GAP.@wrap IsAlmostSimpleGroup(x::Any)::Bool
GAP.@wrap IsAlternatingForm(x::Any)::Bool
GAP.@wrap IsAlternatingGroup(x::Any)::Bool
GAP.@wrap IsAssocWord(x::Any)::Bool
GAP.@wrap IsAtlasCharacterTable(x::GapObj)::Bool
GAP.@wrap IsBiCoset(x::Any)::Bool
GAP.@wrap IsBijective(x::Any)::Bool
GAP.@wrap IsBool(x::Any)::Bool
GAP.@wrap IsCanonicalBasisAlgebraicExtension(x::GapObj)::Bool
GAP.@wrap IsChar(x::Any)::Bool
GAP.@wrap IsCharacteristicSubgroup(x::Any, y::Any)::Bool
GAP.@wrap IsCheapConwayPolynomial(x::Any)::Bool
GAP.@wrap IsCheapConwayPolynomial(x::Any, y::Any)::Bool
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
GAP.@wrap IsDomain(x::Any)::Bool
GAP.@wrap IsDoneIterator(x::Any)::Bool
GAP.@wrap IsDuplicateTable(x::Any)::Bool
GAP.@wrap IsElementaryAbelian(x::Any)::Bool
GAP.@wrap IsElementOfFpGroupFamily(x::GapObj)::Bool
GAP.@wrap IsEmpty(x::Any)::Bool
GAP.@wrap IsFFE(x::Any)::Bool
GAP.@wrap IsFFECollCollColl(x::Any)::Bool
GAP.@wrap IsField(x::Any)::Bool
GAP.@wrap IsFinite(x::Any)::Bool
GAP.@wrap IsFiniteDimensional(x::Any)::Bool
GAP.@wrap IsFinitelyGeneratedGroup(x::Any)::Bool
GAP.@wrap IsFpGroup(x::GapObj)::Bool
GAP.@wrap IsFreeGroup(x::GapObj)::Bool
GAP.@wrap IsGroupOfAutomorphisms(x::Any)::Bool
GAP.@wrap IsHandledByNiceMonomorphism(x::Any)::Bool
GAP.@wrap IsHermitianForm(x::Any)::Bool
GAP.@wrap IsHomogeneousList(x::Any)::Bool
GAP.@wrap IsInjective(x::Any)::Bool
GAP.@wrap IsInnerAutomorphism(x::Any)::Bool
GAP.@wrap IsInt(x::Any)::Bool
GAP.@wrap IsIntegers(x::Any)::Bool
GAP.@wrap IsIrreducibleCharacter(x::Any)::Bool
GAP.@wrap IsLeaf(x::GapObj)::Bool
GAP.@wrap IsLetterAssocWordRep(x::Any)::Bool
GAP.@wrap IsLetterWordsFamily(x::Any)::Bool
GAP.@wrap IsLibTomRep(x::Any)::Bool
GAP.@wrap IsLieAlgebra(x::Any)::Bool
GAP.@wrap IsLieObjectCollection(x::Any)::Bool
GAP.@wrap IsList(x::Any)::Bool
GAP.@wrap IsMatrix(x::GapObj)::Bool
GAP.@wrap IsMatrixGroup(x::GapObj)::Bool
GAP.@wrap IsMatrixObj(x::GapObj)::Bool
GAP.@wrap IsMatrixOrMatrixObj(x::Any)::Bool
GAP.@wrap IsNaturalAlternatingGroup(x::Any)::Bool
GAP.@wrap IsNaturalSymmetricGroup(x::Any)::Bool
GAP.@wrap IsNilpotentGroup(x::Any)::Bool
GAP.@wrap IsNormal(x::Any, y::Any)::Bool
GAP.@wrap IsNumberField(x::Any)::Bool
GAP.@wrap IsNumberFieldByMatrices(x::Any)::Bool
GAP.@wrap IsomorphicSubgroups(x::GapObj, y::GapObj)::GapObj
GAP.@wrap IsomorphismFpGroup(x::GapObj)::GapObj
GAP.@wrap IsomorphismFpGroupByGenerators(x::GapObj, y::GapObj)::GapObj
GAP.@wrap IsomorphismFpGroupByPcgs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap IsOne(x::Any)::Bool
GAP.@wrap IsPcGroup(x::Any)::Bool
GAP.@wrap IsPcpElement(x::Any)::Bool
GAP.@wrap IsPcpGroup(x::Any)::Bool
GAP.@wrap IsPerfectGroup(x::Any)::Bool
GAP.@wrap IsPermGroup(x::Any)::Bool
GAP.@wrap IsPGroup(x::Any)::Bool
GAP.@wrap IsPolynomial(x::Any)::Bool
GAP.@wrap IsPolynomialRing(x::Any)::Bool
GAP.@wrap IsPrimeField(x::Any)::Bool
GAP.@wrap IsPrimitive(x::Any)::Bool
GAP.@wrap IsPrimitive(x::Any, y::Any)::Bool
GAP.@wrap IsQuasisimpleGroup(x::Any)::Bool
GAP.@wrap IsQuaternionGroup(x::Any)::Bool
GAP.@wrap IsRationals(x::Any)::Bool
GAP.@wrap IsReady(x::GapObj)::Bool
GAP.@wrap IsRecogInfoForAlmostSimpleGroup(x::GapObj)::Bool
GAP.@wrap IsRecogInfoForSimpleGroup(x::GapObj)::Bool
GAP.@wrap IsRecord(x::Any)::Bool
GAP.@wrap IsRegular(x::Any)::Bool
GAP.@wrap IsRegular(x::Any, y::Any)::Bool
GAP.@wrap IsSemiRegular(x::Any)::Bool
GAP.@wrap IsSemiRegular(x::Any, y::Any)::Bool
GAP.@wrap IsSet(x::Any)::Bool
GAP.@wrap IsSimpleGroup(x::Any)::Bool
GAP.@wrap IsSingularForm(x::Any)::Bool
GAP.@wrap IsSolvableGroup(x::Any)::Bool
GAP.@wrap IsSporadicSimpleGroup(x::Any)::Bool
GAP.@wrap IsString(x::Any)::Bool
GAP.@wrap IsSubgroupFpGroup(x::Any)::Bool
GAP.@wrap IsSubset(x::Any, y::Any)::Bool
GAP.@wrap IsSupersolvableGroup(x::Any)::Bool
GAP.@wrap IsSurjective(x::Any)::Bool
GAP.@wrap IsSyllableAssocWordRep(x::Any)::Bool
GAP.@wrap IsSyllableWordsFamily(x::Any)::Bool
GAP.@wrap IsSymmetricForm(x::Any)::Bool
GAP.@wrap IsSymmetricGroup(x::Any)::Bool
GAP.@wrap IsTableOfMarks(x::Any)::Bool
GAP.@wrap IsTransitive(x::Any)::Bool
GAP.@wrap IsTransitive(x::Any, y::Any)::Bool
GAP.@wrap IsTrivial(x::Any)::Bool
GAP.@wrap IsUnivariatePolynomialRing(x::Any)::Bool
GAP.@wrap IsWholeFamily(x::Any)::Bool
GAP.@wrap IsZero(x::Any)::Bool
GAP.@wrap IsZmodnZObj(x::Any)::Bool
GAP.@wrap IsZmodnZObjNonprimeCollection(x::Any)::Bool
GAP.@wrap Iterator(x::Any)::GapObj
GAP.@wrap KernelOfCharacter(x::GapObj, y::GapObj)::GapObj
GAP.@wrap KernelRecogNode(x::GapObj)::GapObj
GAP.@wrap LargestMovedPoint(x::Any)::Int
GAP.@wrap LatticeByCyclicExtension(x::GapObj, y::GapObj)::GapObj
GAP.@wrap LeadingExponent(x::GapObj)::GapInt
GAP.@wrap LeadingExponentOfPcElement(x::GapObj, y::GapObj)::GapInt
GAP.@wrap LeftActingDomain(x::GapObj)::GapObj
GAP.@wrap LetterRepAssocWord(x::GapObj)::GapObj
GAP.@wrap LibInfoCharacterTable(x::GapObj)::GapObj
GAP.@wrap LieAlgebraByStructureConstants(x::GapObj, y::GapObj)::GapObj
GAP.@wrap LinearCharacters(x::GapObj)::GapObj
GAP.@wrap LinearCombination(x::GapObj, y::GapObj)::GapObj
GAP.@wrap LinearCombinationPcgs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap ListPerm(x::GapObj)::GapObj
GAP.@wrap MappedWord(x::GapObj, y::GapObj, z::GapObj)::GAP.Obj
GAP.@wrap MarksTom(x::GapObj)::GapObj
GAP.@wrap MatScalarProducts(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap MatTom(x::GapObj)::GapObj
GAP.@wrap Maxes(x::GapObj)::GapObj
GAP.@wrap MinimalGeneratingSet(x::GapObj)::GapObj
GAP.@wrap MinimalPolynomial(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap mod(x::Any, y::Any)::GAP.Obj
GAP.@wrap NameFunction(x::GapObj)::GapObj
GAP.@wrap NamesOfFusionSources(x::GapObj)::GapObj
GAP.@wrap NextIterator(x::GapObj)::Any
GAP.@wrap NiceGens(x::GapObj)::GapObj
GAP.@wrap NormalClosure(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Normalizer(x::GapObj, y::GapObj)::GapObj
GAP.@wrap NormalSubgroupClasses(x::GapObj, y::GAP.Obj)::GapObj
GAP.@wrap NrCols(x::GapObj)::Int
GAP.@wrap NrConjugacyClasses(x::Any)::GapInt
GAP.@wrap NrRows(x::GapObj)::Int
GAP.@wrap NumberColumns(x::GapObj)::Int
GAP.@wrap NumberRows(x::GapObj)::Int
GAP.@wrap NumeratorRat(x::Any)::GapInt
GAP.@wrap ObjByExtRep(x::GapObj, y::GapObj)::GapObj
GAP.@wrap One(x::Any)::GAP.Obj
GAP.@wrap OnIndeterminates(x::GapObj, y::GapObj)::GapObj
GAP.@wrap OnLines(x::GapObj, y::GapObj)::GapObj
GAP.@wrap OnSets(x::GapObj, y::GapObj)::GapObj
GAP.@wrap OnSetsSets(x::GapObj, y::GapObj)::GapObj
GAP.@wrap OnTuples(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Order(x::Any)::GapInt
GAP.@wrap OrthogonalComponents(x::GapObj, y::GapObj, z::GapInt)::GapObj
GAP.@wrap PcElementByExponentsNC(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Pcgs(x::GapObj)::GapObj
GAP.@wrap PCore(x::GapObj, y::GapInt)::GapObj
GAP.@wrap PcpElementByExponentsNC(x::GapObj, y::GapObj)::GapObj
GAP.@wrap PcpGroupByCollectorNC(x::GapObj)::GapObj
GAP.@wrap PermList(x::GapObj)::GapObj
GAP.@wrap PermutationCharacter(x::GapObj, y::GapObj)::GapObj
GAP.@wrap Permuted(x::GapObj, y::GapObj)::GapObj
GAP.@wrap PolynomialByExtRep(x::GapObj, y::GapObj)::GapObj
GAP.@wrap PolynomialRing(x::GapObj)::GapObj
GAP.@wrap PolynomialRing(x::GapObj, y::Int)::GapObj
GAP.@wrap PossibleClassFusions(x::GapObj, y::GapObj)::GapObj
GAP.@wrap PossibleClassFusions(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap POW(x::GAP.Obj, y::GAP.Obj)::GAP.Obj
GAP.@wrap PowerMap(x::GapObj, y::Int)::GapObj
GAP.@wrap PowerMap(x::GapObj, y::Int, z::Int)::Int
GAP.@wrap PrimeBlocks(x::GapObj, y::Int)::GapObj
GAP.@wrap PrimePGroup(x::GapObj)::GapInt
GAP.@wrap PrimitiveElement(x::GapObj)::GapObj
GAP.@wrap Projection(x::GapObj)::GapObj
GAP.@wrap Projection(x::GapObj, i::Int)::GapObj
GAP.@wrap QUO(x::GAP.Obj, y::GAP.Obj)::GAP.Obj
GAP.@wrap Random(x::GapObj, y::GapObj)::GAP.Obj
GAP.@wrap Range(x::GapObj)::GapObj
GAP.@wrap RecognizeGroup(x::GapObj)::GapObj
GAP.@wrap ReduceCoeffs(x::GapObj, y::GapObj)
GAP.@wrap RelativeOrder(x::GapObj)::GapInt
GAP.@wrap RelativeOrderOfPcElement(x::GapObj, y::GapObj)::GapInt
GAP.@wrap RelativeOrders(x::GapObj)::GapObj
GAP.@wrap RelatorsOfFpGroup(x::GapObj)::GapObj
GAP.@wrap Representative(x::GapObj)::GAP.Obj
GAP.@wrap RepresentativeAction(x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap RepresentativeTom(x::GapObj, y::Int)::GapObj
GAP.@wrap RestrictedMapping(x::GapObj, y::GapObj)::GapObj
GAP.@wrap RightCoset(x::GapObj, y::GapObj)::GapObj
GAP.@wrap RightTransversal(x::GapObj, y::GapObj)::GapObj
GAP.@wrap RootSystem(x::GapObj)::GapObj
GAP.@wrap ScalarProduct(x::GapObj, y::GapObj, z::GapObj)::GAP.Obj
GAP.@wrap SchurIndexByCharacter(x::GapObj, y::GapObj, z::GapObj)::GAP.Obj
GAP.@wrap Set(x::GapObj)::GapObj
GAP.@wrap SetConjugacyClasses(x::Any, y::Any)::Nothing
GAP.@wrap SetIsIrreducibleCharacter(x::Any, y::Bool)::Nothing
GAP.@wrap SetMaximalAbelianQuotient(x::Any, y::Any)::Nothing
GAP.@wrap SetSize(x::Any, y::Any)::Nothing
GAP.@wrap SetUnderlyingGroup(x::Any, y::Any)::Nothing
GAP.@wrap ShrinkRowVector(x::GapObj)::Nothing
GAP.@wrap SignPerm(x::GapObj)::Int
GAP.@wrap SignPermGroup(x::GapObj)::Int
GAP.@wrap Size(x::Any)::GapInt
GAP.@wrap SizeOfFieldOfDefinition(x::GapObj, y::GapInt)::GapInt
GAP.@wrap SizesCentralizers(x::GapObj)::GapObj
GAP.@wrap SizesConjugacyClasses(x::GapObj)::GapObj
GAP.@wrap SLPforElement(x::GapObj, y::GapObj)::GAP.Obj
GAP.@wrap SmallestMovedPoint(x::Any)::GapInt
GAP.@wrap Source(x::GapObj)::GapObj
GAP.@wrap Sqrt(x::Int64)::GAP.Obj
GAP.@wrap Stabilizer(v::GapObj, w::Any, x::GapObj, y::GapObj, z::GapObj)::GapObj
GAP.@wrap Stabilizer(x::GapObj, y::Any, z::GapObj)::GapObj
GAP.@wrap StringViewObj(x::Any)::GapObj
GAP.@wrap StructureConstantsTable(x::GapObj)::GapObj
GAP.@wrap StructureDescription(x::GapObj)::GapObj
GAP.@wrap SubgroupNC(x::GapObj, y::GapObj)::GapObj
GAP.@wrap SubsTom(x::GapObj)::GapObj
GAP.@wrap SylowSubgroup(x::GapObj, y::GapInt)::GapObj
GAP.@wrap SymmetricParts(x::GapObj, y::GapObj, z::GapInt)::GapObj
GAP.@wrap Symmetrizations(x::GapObj, y::GapObj, z::GapInt)::GapObj
GAP.@wrap SymplecticComponents(x::GapObj, y::GapObj, z::GapInt)::GapObj
GAP.@wrap TableOfMarks(x::GapObj)::GapObj
GAP.@wrap Trace(x::GapObj, y::GAP.Obj)::GAP.Obj
GAP.@wrap TrivialCharacter(x::GapObj)::GapObj
GAP.@wrap UnderlyingElement(x::GapObj)::GapObj
GAP.@wrap UnderlyingGroup(x::GapObj)::GapObj
GAP.@wrap UnderlyingRingElement(x::GapObj)::GapObj
GAP.@wrap UnivariatePolynomialByCoefficients(x::GapObj, y::GapObj, z::Int)::GapObj
GAP.@wrap Value(x::GapObj, y::Any)::Any
GAP.@wrap Value(x::GapObj, y::Any, z::Any)::Any
GAP.@wrap Value(x::GapObj, y::GapObj, z::GapObj, a::Any)::Any
GAP.@wrap ValueGlobal(x::GapObj)::GAP.Obj
GAP.@wrap ValuesOfClassFunction(x::GapObj)::GapObj
GAP.@wrap WeylGroup(x::GapObj)::GapObj
GAP.@wrap WeylOrbitIterator(x::GapObj, y::GapObj)::GapObj
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
