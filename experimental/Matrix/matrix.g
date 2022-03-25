BindGlobal("HelloWorld", function()
    return "Hello World!";
end);

BindGlobal("JuliaMatrixFamily", NewFamily("JuliaMatrixFamily"));

DeclareRepresentation(
  "IsJuliaMatrixRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsMatrixObj);
  
JuliaMatrixIdentificationType :=
    NewType( JuliaMatrixFamily, IsJuliaMatrixRep);
    
#DeclareOperation( "MakeJuliaMatrixRep", [ IsJuliaObject] );

#InstallMethod( MakeJuliaMatrixRep, "", [IsJuliaObject],
#    function(m)
#        return Objectify(JuliaMatrixIdentificationType, rec(m := m));
#    end
#);

BindGlobal("MakeJuliaMatrixRep",function(m)
  return Objectify(JuliaMatrixIdentificationType,
               rec(m := m));
end);

#DeclareOperation( "\*", [ IsJuliaMatrixRep, IsJuliaMatrixRep ] );
#DeclareOperation( "\/", [ IsJuliaMatrixRep, IsJuliaMatrixRep ] );

#InstallMethod( \*,
#    "for two IsJuliaMatrixRep",
#    [ IsJuliaMatrixRep, IsJuliaMatrixRep ],
#    function( m1, m2 )
#        return MakeJuliaMatrixRep(m1!.m * m2!.m);
#    end );

#InstallMethod( \/,
#    "for two IsJuliaMatrixRep",
#    [ IsJuliaMatrixRep, IsJuliaMatrixRep ],
#    function( m1, m2 )
#        return m1 * Inverse(m2);
#    end );
    
InstallOtherMethod( NumberRows, [IsJuliaMatrixRep], m -> Julia.nrows(m!.m));

InstallOtherMethod( NumberColumns, [IsJuliaMatrixRep], m -> Julia.ncols(m!.m));

InstallOtherMethod( BaseDomain, [IsJuliaMatrixRep], m -> Julia.base_ring(m!.m));

InstallOtherMethod( \[\,\],
    [ IsJuliaMatrixRep, IsPosInt and IsSmallIntRep,
                       IsPosInt and IsSmallIntRep ],
    function( m, i, j )
      return (m!.m)[i,j];
    end );
    
InstallOtherMethod( \*, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return MakeJuliaMatrixRep(m1!.m * m2!.m);
    end );
    
InstallOtherMethod( \/, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return m1 * Inverse(m2);
    end);
    
InstallMethod( InverseMutable, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(InverseOp(m!.m)));

InstallMethod( InverseImmutable, [IsJuliaMatrixRep], m -> MakeImmutable(MakeJuliaMatrixRep(InverseOp(m!.m))));

InstallOtherMethod( AdditiveInverseOp, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(AdditiveInverseOp(m!.m)));

InstallOtherMethod(One, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(OneOp(m!.m)));

InstallOtherMethod(Zero, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(ZeroOp(m!.m)));

InstallOtherMethod(\+, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return MakeJuliaMatrixRep(m1!.m + m2!.m);
    end);
    
InstallOtherMethod(\-, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return MakeJuliaMatrixRep(m1!.m - m2!.m);
    end);

InstallOtherMethod(\^, [IsJuliaMatrixRep, IsInt and IsSmallIntRep], function( m1, k )
        if k < 0 then
            return MakeJuliaMatrixRep( (Inverse(m1)!.m)^(-k) );
        else
            return MakeJuliaMatrixRep( (m1!.m)^k );
        fi;
    end);
    
InstallOtherMethod(\=, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return m1!.m = m2!.m;
    end);
    
InstallOtherMethod(Display, [IsJuliaMatrixRep], 5000, function( m )
    local i,j;
        Print("[ ");
        for i in [1..NumberRows(m)] do
            Print("[ ");
            for j in [1..NumberColumns(m)] do
                Print(m[i,j]);
                if (j < NumberColumns(m)) then
                    Print(", ");
                fi;
            od;
            if (i < NumberRows(m)) then
                Print(" ] \n");
                Print("  ");
            else
                Print(" ] ] \n");
            fi;
        od;
    end);

InstallOtherMethod(Unpack, [IsJuliaMatrixRep], function(m)
    local v, i, j;
    v := EmptyPlist(NumberRows(m)*NumberColumns(m));
    for i in [1..NumberRows(m)] do
        for j in [1..NumberColumns(m)] do
            v[j+NumberColumns(m)*(i-1)] := m[i,j];
        od;
    od;
    return v;
end);

# TODO: Not working. JuliaMatrixRep is immutable
InstallOtherMethod( \[\,\]\:\=,
    [ IsJuliaMatrixRep and IsMutable, IsPosInt and IsSmallIntRep,
                       IsPosInt and IsSmallIntRep, IsObject ],
    function( m, i, j, val )
    local mat;
        mat := m!.m;
        mat[i,j] := val;
    end );

# TODO: Not working
InstallOtherMethod(\<, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return m1!.m < m2!.m;
    end);
