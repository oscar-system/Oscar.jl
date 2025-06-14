BindGlobal( "QQBarFieldElementFam",
  NewFamily( "QQBarFieldElement", IsQQBarFieldElement ) );

SetIsUFDFamily( QQBarFieldElementFam, true );

DeclareRepresentation( "IsQQBarFieldElementRep", IsAttributeStoringRep );

BindGlobal( "QQBarFieldType",
  NewType( QQBarFieldElementFam,
           IsQQBarFieldElement and IsQQBarFieldElementRep ) );

BindGlobal( "QQBarField_Julia", Oscar_jl.algebraic_closure( Oscar_jl.QQ ) );

BindGlobal( "QQBarFieldElementMatrixType",
  Oscar_jl.OscarInterfaceConstants._Matrix_QQBarFieldElem );

InstallMethod( _QQBarFieldElement,
  [ "IsRat" ],
  x -> Oscar_jl.QQBarFieldElem( Oscar_jl.QQ( x ) ) );

InstallMethod( _QQBarFieldElement,
  [ "IsCyc" ],
  function( x )
  local iso;

  iso:= IsoGapOscar( Cyclotomics );
  return QQBarField_Julia( iso( x ) );
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
  x -> QQBarFieldElement( Oscar_jl.inv( JuliaPointer( x ) ) ) );

InstallMethod( ComplexConjugate,
  [ "IsQQBarFieldElement" ],
  x -> QQBarFieldElement( Julia.conj( JuliaPointer( x ) ) ) );

InstallMethod( Sqrt,
  [ "IsQQBarFieldElement" ],
  function( x )
  return QQBarFieldElement( Julia.sqrt( JuliaPointer( x ) ) );
  end );


InstallMethod( DefaultFieldOfMatrix,
  [ "IsMatrix and IsQQBarFieldElementCollColl" ],
  M -> QQBarField );

BindGlobal( "_DisplayStringQQBarMatrixElement", function( x )
  local str, pos;

  str:= String( x );
  if ValueOption( "short" ) = true then
    pos:= Position( str, ' ' );
    str:= str{ [ pos + 1 .. Position( str, ' ', pos ) - 1 ] };
  fi;
  return str;
  end );

InstallMethod( Display,
  [ "IsMatrix and IsQQBarFieldElementCollColl" ],
  function( M )
  local m, n, F, strings, w, z, zstr, row, x;

  m:= NrRows( M );
  n:= NrCols( M );
  if m = 0 or n = 0 then
    TryNextMethod();
  fi;
  F:= DefaultFieldOfMatrix( M );

  Print( m, "x", n, " matrix over ", F, ":\n" );
  strings:= List( M, row -> List( row, _DisplayStringQQBarMatrixElement ) );

  w:= Maximum( List( strings, row -> Maximum( List( row, Length ) ) ) ) + 1;
  z:= _DisplayStringQQBarMatrixElement( Zero( F ) );
  zstr:= String( ".", w );

  for row in strings do
    for x in row do
      if x = z then
        Print( zstr );
      else
        Print( String( x, w ) );
      fi;
    od;
    Print( "\n" );
  od;
  end );

InstallMethod( DefaultFieldByGenerators,
  [ "IsList and IsQQBarFieldElementCollection" ],
  L -> QQBarField );

# The `rand` method for `QQBarFieldElem` in Oscar is based on
# `_flint_rand_states`, we cannot use a GAP random source for it,
# thus `InstallMethodWithRandomSource` cannot be used.
InstallMethod( Random,
  [ "IsQQBarField" ],
  F -> QQBarFieldElement( CallJuliaFunctionWithKeywordArguments( Oscar_jl.rand,
                            [ QQBarField_Julia ],
                            rec( degree:= 4, bits:= 5 ) ) ) );

InstallOtherMethod( AbsoluteValue,
  [ "IsQQBarFieldElement" ], SUM_FLAGS,
  x -> QQBarFieldElement( Julia.abs( JuliaPointer( x ) ) ) );

InstallMethod( Eigenvalues,
  [ "IsQQBarField", "IsMatrix and IsQQBarFieldElementCollColl" ],
  function( F, M )
  local MM;

  MM:= Oscar_jl.matrix( QQBarField_Julia,
         GAPToJulia( QQBarFieldElementMatrixType,
                     List( M, row -> List( row, JuliaPointer ) ) ) );
  return List( JuliaToGAP( IsList, Oscar_jl.eigenvalues( MM ) ),
               QQBarFieldElement );
  end);

BindGlobal( "IsRealQQBarFieldElement",
  x -> IsQQBarFieldElement( x ) and Oscar_jl.is_real( JuliaPointer( x ) ) );
