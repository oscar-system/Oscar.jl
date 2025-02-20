BindGlobal("ReplaceGapFunc", function(name, func)
  local orig_name;
  # Store a copy of the original value, but only if that copy does not yet
  # exist. This ensures we don't overwrite it during a call to `Reread`.
  orig_name := Concatenation("_ORIG_", name);
  if not IsBoundGlobal(orig_name) then
    BindGlobal(orig_name, ValueGlobal(name));
  fi;
  MakeReadWriteGlobal(name);
  UnbindGlobal(name);
  BindGlobal(name, func);
end);

############################################################################

InstallOtherMethod( ImagesRepresentative,
    [ IsActionHomomorphism and HasJuliaData, IsMultiplicativeElementWithInverse ],
function( hom, elm )
    local data;
    data:= JuliaData( hom );
    return Oscar_jl.permutation(data[1], Oscar_jl.group_element(data[2], elm)).X;
end );

InstallMethod( RestrictedMapping,
    CollFamSourceEqFamElms,
    [ IsActionHomomorphism and HasJuliaData, IsGroup ],
function( hom, H )
    local data, OscarG, xset, Omega, Hgens, Hacts, OscarH, res;
    data:= JuliaData( hom ); # the Oscar G-set and the acting Oscar group G
    OscarG:= data[2]; # the acting Oscar group G
    xset:= UnderlyingExternalSet( hom );
    Omega:= HomeEnumerator( xset ); # the set of Oscar objects
    Hgens:= GeneratorsOfGroup( H ); # GAP generators of H
    Hacts:= List( Hgens, x -> Oscar_jl.group_element( OscarG, x ) ); # corresponding Oscar generators of H
    OscarH:= Oscar_jl._as_subgroup_bare( OscarG, H );
    res:= ActionHomomorphism( H, Omega, Hgens, Hacts, FunctionAction( xset ) );
    SetJuliaData( res, [ data[1], OscarH ] );
    return res;
end );

############################################################################

# Install the methods for `IsoGapOscar`,
# in order to make Oscar's `iso_gap_oscar` work.
Perform( Oscar._iso_gap_oscar_methods,
         function( pair )
           InstallMethod( IsoGapOscar, [ JuliaToGAP( IsString, pair[1] ) ],
               Oscar.GAP.WrapJuliaFunc( pair[2] ) );
         end );

############################################################################

# GAP supports classical matrix groups over the ring `R` in the following
# cases.
# - `descr` one of "GL", "SL":
#   finite fields, Z, residue rings of Z,
# - `descr` "Sp":
#   finite fields, residue rings Z/p^n Z for primes p,
# - `descr` one of "GO", "GO+", "GO-", "SO", "SO+", "SO-",
#   "Omega", "Omega+", "Omega-":
#   finite fields, residue rings Z/p^n Z for odd primes p,
# - `descr` one of "GU", "SU":
#   finite fields.
BindGlobal( "IsBaseRingSupportedForClassicalMatrixGroup", function( R, descr )
  if descr in [ "GL", "SL" ] then
    return ( IsField( R ) and IsFinite( R ) ) or
           IsIntegers( R ) or
           ( IsZmodnZObjNonprimeCollection( R ) and
             HasIsWholeFamily( R ) and IsWholeFamily( R ) );
  elif descr = "Sp" then
    return ( IsField( R ) and IsFinite( R ) ) or
           ( IsZmodnZObjNonprimeCollection( R ) and
             HasIsWholeFamily( R ) and IsWholeFamily( R ) and
             IsPrimePowerInt( Size( R ) ) );
  elif descr in [ "GO", "GO+", "GO-", "SO", "SO+", "SO-",
                  "Omega", "Omega+", "Omega-" ] then
    return ( IsField( R ) and IsFinite( R ) ) or
           ( IsZmodnZObjNonprimeCollection( R ) and
             HasIsWholeFamily( R ) and IsWholeFamily( R ) and
             IsOddInt( Size( R ) ) and IsPrimePowerInt( Size( R ) ) );
  elif descr in [ "GU", "SU" ] then
    return IsField( R ) and IsFinite( R );
  fi;
  return false;
end );

############################################################################

BindGlobal("_OSCAR_GroupElem", Oscar.AbstractAlgebra.GroupElem);

# the following code ensures that `GAP.Globals.Group` accepts a GAP list
# of Oscar group elements, as there is a GAP `OrbitStabilizerAlgorithm` method
# which sometimes does that when called by our `stabilize` method.
BindGlobal("JuliaObjectCollectionsFamily", CollectionsFamily(JuliaObjectFamily));
InstallMethod( IsGeneratorsOfMagmaWithInverses,
    "for a list or collection of Julia objects",
    f -> f= JuliaObjectCollectionsFamily,
    [ IsListOrCollection ],
    gens -> IsCollection( gens ) and
       ForAll( gens, x -> Julia.isa(x, _OSCAR_GroupElem) ) );

############################################################################

InstallMethod( GroupGeneratorsDefinePresentation,
   [ "IsPcGroup" ],
   G -> GeneratorsOfGroup(G) = FamilyPcgs(G) );

InstallMethod( GroupGeneratorsDefinePresentation,
   [ "IsPcpGroup" ],
   function( G )
     local n, Ggens, i, w;

     n:= One( G )!.collector![ PC_NUMBER_OF_GENERATORS ];
     Ggens:= GeneratorsOfGroup( G );
     if Length( Ggens ) <> n then
       return false;
     fi;
     for i in [ 1 .. n ] do
       w:= Ggens[i]!.word;
       if not ( Length(w) = 2 and w[1] = i and w[2] = 1 ) then
         return false;
       fi;
     od;
     return true;
   end );

############################################################################


Perform( Oscar._GAP_serializations,
         function( entry )
           InstallMethod( SerializeInOscar,
             [ JuliaToGAP( IsString, entry[1] ), "IsObject" ],
             GAP_jl.WrapJuliaFunc( entry[2] ) );
         end );

Perform( Oscar._GAP_deserializations,
         function( entry )
           InstallMethod( DeserializeInOscar,
             [ JuliaToGAP( IsString, entry[1] ), "IsObject", "IsObject" ],
             GAP_jl.WrapJuliaFunc( entry[2] ) );
         end );

############################################################################

# The following can be removed as soon as the CTblLib package provides it
# (not yet in CTblLib v1.3.9).
if not IsBound( IsAtlasCharacterTable ) then
  DeclareProperty( "IsAtlasCharacterTable", IsNearlyCharacterTable );

  InstallMethod( IsAtlasCharacterTable,
    [ "IsOrdinaryTable" ],
    tbl -> PositionSublist( InfoText( tbl ),
                            "origin: ATLAS of finite groups" ) <> fail );

  InstallMethod( IsAtlasCharacterTable,
    [ "IsBrauerTable" ],
    tbl -> IsAtlasCharacterTable( OrdinaryCharacterTable( tbl ) ) );

  AddSet( CTblLib.SupportedAttributes, "IsAtlasCharacterTable" );

  DatabaseAttributeAddX( CTblLib.Data.IdEnumerator, rec(
    identifier:= "IsAtlasCharacterTable",
    type:= "values",
    name:= "IsAtlasCharacterTable",
    neededAttributes:= [ "InfoText" ],
    create:= function( attr, id )
      local infotext;

      infotext:= attr.idenumerator.attributes.InfoText;
      return PositionSublist( infotext.attributeValue( infotext, id ),
                              "origin: ATLAS of finite groups" ) <> fail;
      end,
  ) );

  CTblLib.ExtendAttributeOfIdEnumeratorExt( "IsAtlasCharacterTable", true );
fi;
