#############################################################################
##
##  ExternalPolymakeCone.gi      JConvex package
##                               Martin Bies
##
##  Copyright 2021               University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake
##
##  Chapter Cones in Polymake
##
#############################################################################


##############################################################################################
##
##  Section GAP category of PolymakeCones
##
##############################################################################################

DeclareRepresentation( "IsPolymakeConeRep", IsPolymakeCone and IsAttributeStoringRep, [ ] );

BindGlobal( "TheFamilyOfPolymakeCones", NewFamily( "TheFamilyOfPolymakeCones" ) );

BindGlobal( "TheTypeOfPolymakeCone", NewType( TheFamilyOfPolymakeCones, IsPolymakeConeRep ) );


##############################################################################################
##
##  Constructors for PolymakeCones
##
##############################################################################################

InstallGlobalFunction( Polymake_ConeByGenerators,
  function( arg )
    local given_rays, i, given_lineality, cone;
    
    if Length( arg ) = 0 then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[1] ) then
        
        given_rays := Filtered( arg[ 1 ], row -> not IsZero( row ) );
        if Length( given_rays ) > 0 then
            
            # non-trivial cone
            return Polymake_ConeByGenerators( arg[ 1 ], [ ] );
            
        else
            
            # received the trivial cone
            cone := rec( rays := [ ],
                         lineality := [ ],
                         number_type := "rational",
                         rep_type := "V-rep" );
            ObjectifyWithAttributes( cone, TheTypeOfPolymakeCone );
            return cone;
            
        fi;
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if ( not IsEmpty( arg[ 1 ] ) ) and not ( IsMatrix( arg[ 1 ] ) ) then
            Error( "Wronge input: The first argument should be a Gap matrix!" );
        fi;
        
        if ( not IsEmpty( arg[ 2 ] ) ) and not ( IsMatrix( arg[ 2 ] ) ) then
            Error( "Wronge input: The second argument should be a Gap matrix!" );
        fi;
        
        given_rays := Filtered( arg[ 1 ], row -> not IsZero( row ) );
        given_lineality := Filtered( arg[ 2 ], row -> not IsZero( row ) );
        cone := rec( rays := given_rays,
                     lineality := given_lineality,
                     number_type := "rational",
                     rep_type := "V-rep" );
        ObjectifyWithAttributes( cone, TheTypeOfPolymakeCone );
        
        return Polymake_CanonicalConeByGenerators( cone );
        
    fi;
    
end );

InstallGlobalFunction( Polymake_ConeFromInequalities,
  function( arg )
    local given_ineqs, given_eqs, cone;
    
    if Length( arg ) = 0 or ForAll( arg, IsEmpty ) then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[ 1 ] ) then
        
        return Polymake_ConeFromInequalities( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if ( not IsEmpty( arg[ 1 ] ) ) and not ( IsMatrix( arg[ 1 ] ) ) then
            Error( "Wronge input: The first argument should be a Gap matrix!" );
        fi;
        
        if ( not IsEmpty( arg[ 2 ] ) ) and not ( IsMatrix( arg[ 2 ] ) ) then
            Error( "Wronge input: The second argument should be a Gap matrix!" );
        fi;
        
        given_ineqs := Filtered( arg[ 1 ], row -> not IsZero( row ) );
        given_eqs := Filtered( arg[ 2 ], row -> not IsZero( row ) );
        cone := rec( inequalities := given_ineqs,
                     equalities := given_eqs,
                     number_type := "rational",
                     rep_type := "H-rep" );
        ObjectifyWithAttributes( cone, TheTypeOfPolymakeCone );
        
        return Polymake_CanonicalConeFromInequalities( cone );
        
    fi;
    
end );


##############################################################################################
##
##  Canonicalize cones
##
##############################################################################################

InstallMethod( Polymake_CanonicalConeByGenerators,
               [ IsPolymakeCone ],
  function( cone )
    local command_string, s, res_string, rays, scaled_rays, i, scale, lineality, scaled_lineality, new_cone;
    
    if cone!.rep_type = "H-rep" then
        
        return fail;
        
    else
        
        # compute rays
        command_string := Concatenation( Polymake_V_Rep_command_string( cone ), ".RAYS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        rays := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational rays - we turn them into integral vectors
        scaled_rays := [];
        for i in [ 1 .. Length( rays ) ] do
            scale := Lcm( List( rays[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_rays, [ scale * rays[ i ] ] );
        od;
        
        # extract lineality
        command_string := Concatenation( Polymake_V_Rep_command_string( cone ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational lineality - we turn them into integral vectors
        scaled_lineality := [];
        for i in [ 1 .. Length( lineality ) ] do
            scale := Lcm( List( lineality[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_lineality, [ scale * lineality[ i ] ] );
        od;
        
        # construct the new cone
        new_cone := rec( rays := scaled_rays,
                         lineality := scaled_lineality,
                         number_type := "rational",
                         rep_type := "V-rep" );
        ObjectifyWithAttributes( new_cone, TheTypeOfPolymakeCone );
        return new_cone;
        
    fi;
    
end );

InstallMethod( Polymake_CanonicalConeFromInequalities,
               [ IsPolymakeCone ],
  function( cone )
    local command_string, s, res_string, ineqs, scaled_ineqs, i, scale, eqs, scaled_eqs, new_cone;
    
    if cone!.rep_type = "V-rep" then
        
        return fail;
        
    else
        
        # compute facets
        command_string := Concatenation( Polymake_H_Rep_command_string( cone ), ".FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational facets - we turn them into integral vectors
        scaled_ineqs := [];
        for i in [ 1 .. Length( ineqs ) ] do
            scale := Lcm( List( ineqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_ineqs, [ scale * ineqs[ i ] ] );
        od;
        
        # compute linear span
        command_string := Concatenation( Polymake_H_Rep_command_string( cone ), ".LINEAR_SPAN" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational linear spans - we turn them into integral vectors
        scaled_eqs := [];
        for i in [ 1 .. Length( eqs ) ] do
            scale := Lcm( List( eqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_eqs, [ scale * eqs[ i ] ] );
        od;
        
        # construct the new cone
        new_cone := rec( inequalities := scaled_ineqs,
                         equalities := scaled_eqs,
                         number_type := "rational",
                         rep_type := "H-rep" );
        ObjectifyWithAttributes( new_cone, TheTypeOfPolymakeCone );
        return new_cone;
        
    fi;
    
end );


##############################################################################################
##
##  Conversion of cones
##
##############################################################################################

InstallMethod( Polymake_V_Rep,
               [ IsPolymakeCone ],
  function( cone )
    local command_string, s, res_string, rays, scaled_rays, i, scale, lineality, scaled_lineality, new_cone;
    
    if cone!.rep_type = "V-rep" then
        return cone;
    else
        
        # compute rays
        command_string := Concatenation( Polymake_H_Rep_command_string( cone ), ".RAYS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        rays := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational rays - we turn them into integral vectors
        scaled_rays := [];
        for i in [ 1 .. Length( rays ) ] do
            scale := Lcm( List( rays[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_rays, [ scale * rays[ i ] ] );
        od;
        
        # extract lineality
        command_string := Concatenation( Polymake_H_Rep_command_string( cone ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational lineality - we turn them into integral vectors
        scaled_lineality := [];
        for i in [ 1 .. Length( lineality ) ] do
            scale := Lcm( List( lineality[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_lineality, [ scale * lineality[ i ] ] );
        od;
        
        # construct the new cone
        new_cone := rec( rays := scaled_rays,
                         lineality := scaled_lineality,
                         number_type := "rational",
                         rep_type := "V-rep" );
        ObjectifyWithAttributes( new_cone, TheTypeOfPolymakeCone );
        return new_cone;
        
    fi;
    
end );

InstallMethod( Polymake_H_Rep,
               [ IsPolymakeCone ],
  function( cone )
    local command_string, s, res_string, ineqs, scaled_ineqs, i, scale, eqs, scaled_eqs, new_cone;
    
    if cone!.rep_type = "H-rep" then
        
        return cone;
        
    else
        
        if cone!.rep_type = "V-rep" and cone!.rays = [] then
            return Polymake_ConeFromInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # compute facets
        command_string := Concatenation( Polymake_V_Rep_command_string( cone ), ".FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational facets - we turn them into integral vectors
        scaled_ineqs := [];
        for i in [ 1 .. Length( ineqs ) ] do
            scale := Lcm( List( ineqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_ineqs, [ scale * ineqs[ i ] ] );
        od;
        
        # compute linear span
        command_string := Concatenation( Polymake_V_Rep_command_string( cone ), ".LINEAR_SPAN" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational linear spans - we turn them into integral vectors
        scaled_eqs := [];
        for i in [ 1 .. Length( eqs ) ] do
            scale := Lcm( List( eqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_eqs, [ scale * eqs[ i ] ] );
        od;
        
        # construct the new cone
        new_cone := rec( inequalities := scaled_ineqs,
                         equalities := scaled_eqs,
                         number_type := "rational",
                         rep_type := "H-rep" );
        ObjectifyWithAttributes( new_cone, TheTypeOfPolymakeCone );
        return new_cone;
        
    fi;
    
end );


##############################################################################################
##
##  Attributes of cones
##
##############################################################################################

InstallMethod( Polymake_AmbientSpaceDimension,
              "finding the dimension of the ambient space of the cone",
              [ IsPolymakeCone ],
  function( cone )
    
    return Length( Polymake_V_Rep( cone )!.rays[1] );
    
end );


InstallMethod( Polymake_Dimension,
              " returns the dimension of the cone",
            [ IsPolymakeCone ],
  function( cone )
    local command_string, s;
    
    if Polymake_IsEmpty( cone ) then 
        return -1;
    fi;
    
    command_string := Concatenation( Polymake_V_Rep_command_string( Polymake_V_Rep( cone ) ), ".CONE_DIM" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


InstallMethod( Polymake_Rays,
              " return the list of generating vertices",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( Polymake_V_Rep( cone )!.rays );
    
end );


InstallMethod( Polymake_Lineality,
              " return the list of generating vertices",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( Polymake_V_Rep( cone )!.lineality );
    
end );


InstallMethod( Polymake_Equalities,
              " return the list of equalities of a cone",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( ( Polymake_H_Rep( cone ) )!.equalities );
    
end );


InstallMethod( Polymake_Inequalities,
              " return the list of inequalities of a cone",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( ( Polymake_H_Rep( cone ) )!.inequalities );
    
end );


InstallMethod( Polymake_RaysInFacets,
              " returns the incident matrix of the rays in the facets",
            [ IsPolymakeCone ],
  function( cone )
    local help_cone, command_string, s, res_string, number_rays, ray_list, i, dummy, j, helper;
    
    # produce command string
    help_cone := Polymake_V_Rep( cone );
    number_rays := Length( help_cone!.rays );
    command_string := Concatenation( Polymake_V_Rep_command_string( help_cone ), ".RAYS_IN_FACETS" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
    
    # process the string
    res_string := SplitString( s, '\n' );
    res_string := List( [ 2 .. Length( res_string ) ], i -> EvalString( res_string[ i ] ) );
    
    # now construct list of rays -- each facets is a list number_rays entries: 0 indicates that the ray generator is not in the facet and a 1 that it is part of the facet
    ray_list := [];
    for i in [ 1 .. Length( res_string ) ] do
        dummy := List( [ 1 .. number_rays ], i -> 0 );
        for j in [ 1 .. Length( res_string[ i ] ) ] do
            dummy[ res_string[ i ][ j ] ] := 1;
        od;
        Append( ray_list, [ dummy ] );
    od;
    
    # return the result
    return ray_list;
    
end );

InstallMethod( Polymake_RaysInFaces,
              " returns the incident matrix of the rays in the faces",
            [ IsPolymakeCone ],
  function( cone )
    local dim, help_cone, rays, rays_in_faces, i, generators, additional_rays, j, converted_additional_rays, pos, pos2, k, help_list;
    
    # check degenerate case
    dim := Polymake_Dimension( cone );
    if ( dim = 2 ) then
        
        rays_in_faces := Polymake_RaysInFacets( cone );
        
    else
        
        # compute a V-representation of the cone
        help_cone := Polymake_V_Rep( cone );
        
        # read-off the rays and the list of rays in the facets
        rays := help_cone!.rays;
        rays_in_faces := Polymake_RaysInFacets( help_cone );
        
        # construct the facets as cones, compute their facets and add them
        # this is thus a recursive function, with (presumably) low performance
        for i in [ 1 .. Length( rays_in_faces ) ] do
            
            # construct the i-th facet and compute the rays in its facets
            generators := List( Filtered( [ 1 .. Length( rays ) ], j -> not IsZero( rays_in_faces[ i ][ j ] ) ), k -> rays[ k ] );
            additional_rays := RaysInFaces( Cone( generators ) );
            
            # convert these lists into a list with the rays above
            converted_additional_rays := [];
            for j in [ 1 .. Length( additional_rays ) ] do
                pos := Positions( additional_rays[ j ], 1 );
                pos2 := List( [ 1 .. Length( pos ) ], k -> Position( rays, generators[ pos[ k ] ] ) );
                help_list := List( [ 1 .. Length( rays ) ], k -> 0 );
                for k in pos2 do
                    help_list[ k ] := 1;
                od;
                converted_additional_rays := Concatenation( converted_additional_rays, [ help_list ] );
            od;
            
            # add this list of rays in the faces of the i-th facet to the above-computed list of rays_in_facets
            rays_in_faces := Concatenation( rays_in_faces, converted_additional_rays  );
            
        od;
        
    fi;
    
    # all all rays, since the overall cone is also a face
    rays_in_faces := Concatenation( rays_in_faces, [ List( [ 1 .. Length( Polymake_Rays( cone ) ) ], i -> 1 ) ] );
    
    # return the result
    return DuplicateFreeList( rays_in_faces );
    
end );


##############################################################################################
##
##  Properties of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_IsEmpty,
               "finding if the cone empty is or not",
               [ IsPolymakeCone ],
  function( cone )
    
    return Length( Polymake_V_Rep( cone )!.rays ) = 0;
    
end );


InstallMethod( Polymake_IsPointed,
               "finding if the cone is pointed or not",
               [ IsPolymakeCone ],
  function( cone )
    local help_cone, rays, lin, command_string, s;
    
    # compute V-representation
    help_cone := Polymake_V_Rep( cone );
    
    # parse the rays into format recognized by Polymake
    command_string := Concatenation( Polymake_V_Rep_command_string( help_cone ), ".POINTED" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
    
    # return the result
    return EvalString( s );
    
end );


##############################################################################################
##
##  Operations with cones
##
##############################################################################################

InstallMethod( Polymake_Intersection,
               "construct intersection of two cones",
               [ IsPolymakeCone, IsPolymakeCone ],
  function( cone1, cone2 )
    local cone1_h, cone2_h, new_ineqs, new_equ;
    
    cone1_h := Polymake_H_Rep( cone1 );
    cone2_h := Polymake_H_Rep( cone2 );
    
    new_ineqs := Concatenation( cone1_h!.inequalities, cone2_h!.inequalities );
    new_equ := Concatenation( cone1_h!.equalities, cone2_h!.equalities );
    
    return Polymake_ConeFromInequalities( new_ineqs, new_equ );
    
end );


##############################################################################################
##
##  Command strings
##
##############################################################################################

InstallMethod( Polymake_V_Rep_command_string,
               "construct command string for V-Representation of Cone in Julia",
               [ IsPolymakeCone ],
  function( cone )
    local rays, lin, command_string;
        
        # check if the given cone is a V-rep
        if not ( cone!.rep_type = "V-rep" ) then
            return "fail";
        fi;
        
        # extract data
        rays := cone!.rays;
        lin := cone!.lineality;
        
        # check for degenerate case
        if Length( rays ) = 0 then
            rays := lin;
        fi;
        
        # prepare rays (always non-zero) for command string
        rays := List( [ 1 .. Length( rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "ConeByGAP4PackageConvex", " = Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( rays, "; " ), "] " );
        
        # if we also have lineality, add them
        lin := cone!.lineality;
        if ( Length( lin ) > 0 ) then
            lin := List( [ 1 .. Length( lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, ", INPUT_LINEALITY = [ ", JoinStringsWithSeparator( lin, "; " ), " ] " );
        fi;
        
        # return command string
        return Concatenation( command_string, ")" );
    
end );

InstallMethod( Polymake_H_Rep_command_string,
               "construct command string for H-Representation of Cone in Julia",
               [ IsPolymakeCone ],
  function( cone )
    local ineqs, eqs, command_string;
    
        # check if the given cone is a V-rep
        if not ( cone!.rep_type = "H-rep" ) then
            return fail;
        fi;
        
        # extract data
        ineqs := cone!.inequalities;
        eqs := cone!.equalities;
        
        # check for degenerate case
        if Length( ineqs ) = 0 then
            ineqs := eqs;
        fi;
        
        # prepare inequalities (always non-zero) for command string
        ineqs := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "ConeByGAP4PackageConvex", " = Polymake.polytope.Cone( INEQUALITIES = [ ", JoinStringsWithSeparator( ineqs, "; " ), " ] " );
        
        # if we have non-zero equalities, add them
        eqs := cone!.equalities;
        if ( Length( eqs ) > 0 ) then
            eqs := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, ", EQUATIONS = [ ", JoinStringsWithSeparator( eqs, "; " ), " ] " );
        fi;
        
        # return command string
        return Concatenation( command_string, ")" );
        
end );
