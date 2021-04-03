#############################################################################
##
##  Tools.gi            Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
## Chapter: Tools
##
#############################################################################

InstallMethod( GeneratingVerticesAndGeneratingRays,
               [ IsList, IsList ],
  function( matrix, linearity )
    local generating_vertices, generating_rays, current, temp, l, i;
    
    generating_vertices:= [  ];
    
    generating_rays:= [  ];
    
    temp := StructuralCopy( matrix );
    
    l := Length( temp );
    
    for i in [ 1..l ] do
      
      current := temp[ i ];
      
      if current[ 1 ] = 1 then
        
        Remove( current, 1 );
        
        if i in linearity then
          
          Add( generating_vertices, current );
          
          Add( generating_vertices, -current );
        
        else
          
          Add( generating_vertices, current );
        
        fi;
       
      else
      
        Remove( current, 1 ); 
        
        if not IsZero( current ) then
          
          if i in linearity then
            
            Add( generating_rays, current );
            
            Add( generating_rays, -current );
          
          else
            
            Add( generating_rays, current );
          
          fi;
        
        else 
         
         Add( generating_vertices, current );
         
        fi;
      
      fi;
    
    od;
    
    return [ generating_vertices, generating_rays ];
    
end );


InstallMethod( InequalitiesAndEqualities,
               [ IsList, IsList ],
  function( matrix, linearity )
    local equalities, inequalities , current, temp, l, i;
    
    inequalities:= [];
    
    equalities:= [];
    
    l:= Length( matrix );
    
    for i in [ 1..l ] do
      
      current:= matrix[ i ];
      
      if i in linearity then
        
        Add( equalities, current );
      
      else
        
        Add( inequalities, current );
       
      fi;
    
    od;
    
  return [ inequalities, equalities ];

end );
