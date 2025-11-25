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
Perform( Oscar_jl._iso_gap_oscar_methods,
         function( pair )
           InstallMethod( IsoGapOscar, [ JuliaToGAP( IsString, pair[1] ) ],
               GAP_jl.WrapJuliaFunc( pair[2] ) );
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

BindGlobal("_OSCAR_GroupElem", Oscar_jl.AbstractAlgebra.GroupElem);

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


Perform( Oscar_jl.Serialization._GAP_type_params,
         function( entry )
           InstallMethod( SerializationInOscarDependentObjects,
             [ JuliaToGAP( IsString, entry[1] ) ],
             GAP_jl.WrapJuliaFunc( entry[2] ) );
         end );

Perform( Oscar_jl.Serialization._GAP_serializations,
         function( entry )
           InstallMethod( SerializeInOscar,
             [ JuliaToGAP( IsString, entry[1] ), "IsObject" ],
             GAP_jl.WrapJuliaFunc( entry[2] ) );
         end );

Perform( Oscar_jl.Serialization._GAP_deserializations,
         function( entry )
           if entry[3] then
             InstallMethod( DeserializeInOscar,
               [ JuliaToGAP( IsString, entry[1] ), "IsObject", "IsObject", "IsObject" ],
               GAP_jl.WrapJuliaFunc( entry[2] ) );
           else
             InstallMethod( DeserializeInOscar,
               [ JuliaToGAP( IsString, entry[1] ), "IsObject", "IsObject" ],
               GAP_jl.WrapJuliaFunc( entry[2] ) );
           fi;
         end );

############################################################################

# Oscar provides the reduction of a finite matrix group over a number field
# to one over a finite field,
# hence these groups can be handled via nice monomorphisms.
#
# In the construction of the `GapObj` of a matrix group over a number field,
# we store the Oscar group in the attribute `JuliaData` and set the flag
# `MayBeHandledByNiceMonomorphism`.
# GAP operations that have a special method for groups in the filter
# `IsHandledByNiceMonomorphism` have also a method for groups in
# `MayBeHandledByNiceMonomorphism`.
# The latter method calls `IsHandledByNiceMonomorphism`,
# and we provide the following method for matrix groups that are in
# `MayBeHandledByNiceMonomorphism` and that store an Oscar group.

InstallMethod( IsHandledByNiceMonomorphism,
    [ "IsMatrixGroup and MayBeHandledByNiceMonomorphism and HasJuliaData" ],
    RankFilter( IsCyclotomicCollCollColl ), # above `IsCyclotomicMatrixGroup`
    function( G )
    if Oscar.is_finite( JuliaData( G ) ) then
      # The `is_finite` call already triggers the computation
      # of the reduction, without error in the case of an infinite group.
      # Make sure that the reduction really gets computed.
      Oscar.isomorphic_group_over_finite_field( JuliaData( G ) );
      Assert( 0, HasNiceMonomorphism( G ) );
      return true;
    fi;
    return false;
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
