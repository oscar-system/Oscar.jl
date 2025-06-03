BindGlobal( "QQBarFieldElementFam",
  NewFamily( "QQBarFieldElement", IsQQBarFieldElement ) );

SetIsUFDFamily( QQBarFieldElementFam, true );

DeclareRepresentation( "IsQQBarFieldElementRep",
  IsJuliaWrapper and IsAttributeStoringRep );

BindGlobal( "QQBarFieldType",
  NewType( QQBarFieldElementFam,
           IsQQBarFieldElement and IsQQBarFieldElementRep ) );

BindGlobal( "QQBarField_Julia", Oscar.algebraic_closure( Oscar.QQ ) );

BindGlobal( "QQBarFieldElementMatrixType",
  Oscar_jl.eval( Julia.Meta.parse( Julia.String( "Matrix{QQBarFieldElem}" ) ) ) );

InstallMethod( _QQBarFieldElement,
  [ "IsRat" ],
  x -> Oscar.QQBarFieldElem( Oscar.QQ( x ) ) );

InstallMethod( _QQBarFieldElement,
  [ "IsCyc" ],
  function( x )
  local iso;

  iso:= Oscar.iso_gap_oscar( Cyclotomics );
  return Oscar.QQBarFieldElem( QQBarField_Julia( iso( x ) ) );
  end );

InstallMethod( _QQBarFieldElement,
  [ "IsJuliaObject" ],
  IdFunc );

BindGlobal( "QQBarFieldElement", function( elm )
  return ObjectifyWithAttributes( rec(), QQBarFieldType,
           JuliaPointer, _QQBarFieldElement( elm ) );
  end );

SetOne( QQBarFieldElementFam, QQBarFieldElement( 1 ) );
SetZero( QQBarFieldElementFam, QQBarFieldElement( 0 ) );


BindGlobal( "QQBarField",
  ObjectifyWithAttributes( rec(),
    NewType( CollectionsFamily( QQBarFieldElementFam ),
             IsAttributeStoringRep and
             IsQQBarField ),
    JuliaPointer, QQBarField_Julia ) );

SetName( QQBarField, "QQBarField" );
SetIsLeftActedOnByDivisionRing( QQBarField, true );
SetSize( QQBarField, infinity );
SetIsFiniteDimensional( QQBarField, false );
SetLeftActingDomain( QQBarField, Rationals );
SetCharacteristic( QQBarField, 0 );
SetPrimeField( QQBarField, Rationals );

SetZero( QQBarField, Zero( QQBarFieldElementFam ) );
SetOne( QQBarField, One( QQBarFieldElementFam ) );


InstallMethod( ViewString,
  [ "IsQQBarFieldElement" ],
  x -> Concatenation( "<", PrintString( JuliaPointer( x  ) ), ">" ) );

InstallMethod( PrintString,
  [ "IsQQBarFieldElement" ],
  x -> Concatenation( "<", PrintString( JuliaPointer( x  ) ), ">" ) );

InstallMethod( String,
  [ "IsQQBarFieldElement" ],
  x -> Concatenation( "<", PrintString( JuliaPointer( x  ) ), ">" ) );

InstallMethod( String,
  [ "IsQQBarField" ],
  x -> "QQBarField" );


InstallMethod( \=,
  [ "IsQQBarFieldElement", "IsQQBarFieldElement" ],
  { x, y } -> JuliaPointer( x ) = JuliaPointer( y ) );

InstallMethod( \=,
  [ "IsQQBarFieldElement", "IsCyc" ],
  { x, y } -> x = QQBarFieldElement( y ) );

InstallMethod( \=,
  [ "IsCyc", "IsQQBarFieldElement" ],
  { x, y } -> QQBarFieldElement( x ) = y );


InstallMethod( \<,
   "for QQBarField elements",
 [ "IsQQBarFieldElement", "IsQQBarFieldElement" ],
 { x, y } -> JuliaPointer( x ) < JuliaPointer( y ) );

InstallMethod( \<,
 [ "IsQQBarFieldElement", "IsCyc" ],
 { x, y } -> x < QQBarFieldElement( y ) );

InstallMethod( \<,
 [ "IsCyc", "IsQQBarFieldElement" ],
 { x, y } -> QQBarFieldElement( x ) < y );

InstallMethod( \in,
  [ "IsObject", "IsQQBarField" ], SUM_FLAGS,
  { x, F } -> IsQQBarFieldElement( x ) );


InstallMethod( \+,
  [ "IsQQBarFieldElement", "IsQQBarFieldElement" ],
  { x, y } -> QQBarFieldElement( JuliaPointer( x ) + JuliaPointer( y ) ) );

InstallMethod( \+,
  [ "IsCyclotomic", "IsQQBarFieldElement" ],
  { x, y } -> QQBarFieldElement( x ) + y );

InstallMethod( \+,
  [ "IsQQBarFieldElement", "IsCyclotomic" ],
  { x, y } -> x + QQBarFieldElement( y ) );


InstallMethod( AdditiveInverseSameMutability,
  [ "IsQQBarFieldElement" ],
  x -> QQBarFieldElement( - JuliaPointer( x ) ) );

InstallMethod( AdditiveInverseMutable,
  [ "IsQQBarFieldElement" ],
  x -> QQBarFieldElement( - JuliaPointer( x ) ) );

InstallMethod( \*,
  [ "IsQQBarFieldElement", "IsQQBarFieldElement" ],
  { x, y } -> QQBarFieldElement( JuliaPointer( x ) * JuliaPointer( y ) ) );

InstallMethod( \*,
  [ "IsQQBarFieldElement", "IsCyc" ],
  { x, y } -> x * QQBarFieldElement( y ) );

InstallMethod( \*,
  [ "IsCyc", "IsQQBarFieldElement" ],
  { x, y } -> QQBarFieldElement( x ) * y );

InstallMethod( InverseOp,
  [ "IsQQBarFieldElement" ],
  x -> QQBarFieldElement( Oscar.inv( x ) ) );

InstallMethod( ComplexConjugate,
  [ "IsQQBarFieldElement" ],
  x -> QQBarFieldElement( Julia.conj( x ) ) );

InstallMethod( Sqrt,
  [ "IsQQBarFieldElement" ],
  function( x )
  return QQBarFieldElement( Julia.sqrt( x ) );
  end );

InstallMethod( DefaultFieldOfMatrix,
  [ "IsMatrix and IsQQBarFieldElementCollColl" ],
  M -> QQBarField );

InstallMethod( DefaultFieldByGenerators,
  [ "IsList and IsQQBarFieldElementCollection" ],
  L -> QQBarField );

# The `rand` method for `QQBarFieldElem` in Oscar is based on
# `_flint_rand_states`, we cannot use a GAP random source for it,
# thus `InstallMethodWithRandomSource` cannot be used.
InstallMethod( Random,
  [ "IsQQBarField" ],
  F -> QQBarFieldElement( CallJuliaFunctionWithKeywordArguments( Oscar.rand,
                            [ QQBarField_Julia ],
                            rec( degree:= 4, bits:= 5 ) ) ) );

InstallOtherMethod( AbsoluteValue,
  [ "IsQQBarFieldElement" ], SUM_FLAGS,
  x -> QQBarFieldElement( Julia.abs( x ) ) );

InstallMethod( Eigenvalues,
  [ "IsQQBarField", "IsMatrix and IsQQBarFieldElementCollColl" ],
  function( F, M )
  local MM;

  MM:= Oscar.matrix( QQBarField_Julia,
         GAPToJulia( QQBarFieldElementMatrixType, M ) );
  return List( JuliaToGAP( IsList, Oscar.eigenvalues( MM ) ),
               QQBarFieldElement );
  end);

BindGlobal( "IsRealQQBarFieldElement",
  x -> IsQQBarFieldElement( x ) and Oscar.is_real( x ) );
