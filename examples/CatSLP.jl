using HomalgProject

import Oscar
using Oscar.StraightLinePrograms

LoadPackage( "RingsForHomalg" )
LoadPackage( "LinearAlgebraForCAP" )
LoadPackage( "GeneralizedMorphismsForCAP" )
SwitchGeneralizedMorphismStandard( g"cospan" )
LoadPackage( "LazyCategories" )


const sentinel = 0x61a5355768bfa46bc93b1f999b51ba9a

# convert the computation represented in x to a `Lazy` SLP
function tolazy(x, dict=IdDict{Any,Lazy}())
    r = get(dict, x, sentinel)
    r !== sentinel && return r

    r =
        if IsCapCategoryObject(x) || IsCapCategoryMorphism(x)
            if HasGenesisOfCellOperation(x)
                @assert HasGenesisOfCellArguments(x)
                operation = Symbol(GenesisOfCellOperation(x))
                operation = getfield(@__MODULE__, operation)
                arguments = GenesisOfCellArguments(x)
                argumentslazy = (tolazy(arguments[i], dict) for i in 1:length(arguments))
                call(argumentslazy...) do args...
                    operation(map(arg -> arg isa Vector ? GapObj(arg) : arg, args)...)
                end
            else
                Lazy(x)
            end
        elseif IsList(x)
            xs = [tolazy(x[i], dict) for i = 1:length(x)]
            list(Lazy, xs)
        elseif x isa Int
            Lazy(x)
        else
            error("invalid")
        end
    dict[x] = r
end

function build_morphism()
    ℚ = HomalgFieldOfRationals( )
    id = HomalgIdentityMatrix( 8, ℚ )
    a = CertainColumns( CertainRows( id, [ 1, 2, 3 ] ), [ 2, 3, 4, 5 ] )
    b = CertainColumns( CertainRows( id, [ 1, 2, 3, 4 ] ), [ 2, 3, 4, 5, 6 ] )
    c = CertainColumns( CertainRows( id, [ 1, 2, 3, 4, 5 ] ), [ 3, 4, 5, 6, 7, 8 ] )
    IsZero( a * b )
    IsZero( b * c )
    IsZero( a * b * c )
    ℚmat = MatrixCategory( ℚ )
    Lazy = LazyCategory( ℚmat, show_evaluation = true )
    a = a / Lazy
    SetLabel( a, "a" )
    b = b / Lazy
    SetLabel( b, "b" )
    c = c / Lazy
    SetLabel( c, "c" )
    d = CokernelProjection( a )
    e = CokernelColift( a, PreCompose( b, c ) )
    f = KernelEmbedding( e )
    g = KernelEmbedding( c )
    h = KernelLift( c, PreCompose( a, b ) )
    i = CokernelProjection( h )
    ff = AsGeneralizedMorphism( f )
    dd = AsGeneralizedMorphism( d )
    bb = AsGeneralizedMorphism( b )
    gg = AsGeneralizedMorphism( g )
    ii = AsGeneralizedMorphism( i )
    ss = PreCompose( [ ff, PseudoInverse( dd ), bb, PseudoInverse( gg ), ii ] )
    HonestRepresentative( ss )
end

function example()
    s = build_morphism()
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        VisualizeInJulia( s )
    else
        Visualize(s)
    end

    l = tolazy(s)
    evaluate(l, [])
end
