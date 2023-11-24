using Documenter

#
# This module only exists to "host" a doctest used by the test suite.
#
module AuxDocTest_DynkinDiagram
@doc raw"""
```jldoctest show_dynkin_diagram.test
julia> using Oscar

julia> show_dynkin_diagram(:A, 1)
1

julia> show_dynkin_diagram(:A, 2)
1 - 2

julia> show_dynkin_diagram(:A, 3)
1 - 2 - 3

julia> show_dynkin_diagram(:A, 4)
1 - 2 - 3 - 4

julia> show_dynkin_diagram(:A, 5)
1 - 2 - 3 - 4 - 5

julia> show_dynkin_diagram(:A, 10)
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 - 10

julia> show_dynkin_diagram(:A, 123)
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 - 10 - 11 - 12 - 13 - 14 - 15 - 16 - 17 - 18 - 19 - 20 - 21 - 22 - 23 - 24 - 25 - 26 - 27 - 28 - 29 - 30 - 31 - 32 - 33 - 34 - 35 - 36 - 37 - 38 - 39 - 40 - 41 - 42 - 43 - 44 - 45 - 46 - 47 - 48 - 49 - 50 - 51 - 52 - 53 - 54 - 55 - 56 - 57 - 58 - 59 - 60 - 61 - 62 - 63 - 64 - 65 - 66 - 67 - 68 - 69 - 70 - 71 - 72 - 73 - 74 - 75 - 76 - 77 - 78 - 79 - 80 - 81 - 82 - 83 - 84 - 85 - 86 - 87 - 88 - 89 - 90 - 91 - 92 - 93 - 94 - 95 - 96 - 97 - 98 - 99 - 100 - 101 - 102 - 103 - 104 - 105 - 106 - 107 - 108 - 109 - 110 - 111 - 112 - 113 - 114 - 115 - 116 - 117 - 118 - 119 - 120 - 121 - 122 - 123

julia> show_dynkin_diagram(:B, 2)
1 >=> 2

julia> show_dynkin_diagram(:B, 3)
1 - 2 >=> 3

julia> show_dynkin_diagram(:B, 4)
1 - 2 - 3 >=> 4

julia> show_dynkin_diagram(:B, 5)
1 - 2 - 3 - 4 >=> 5

julia> show_dynkin_diagram(:B, 10)
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 >=> 10

julia> show_dynkin_diagram(:B, 123)
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 - 10 - 11 - 12 - 13 - 14 - 15 - 16 - 17 - 18 - 19 - 20 - 21 - 22 - 23 - 24 - 25 - 26 - 27 - 28 - 29 - 30 - 31 - 32 - 33 - 34 - 35 - 36 - 37 - 38 - 39 - 40 - 41 - 42 - 43 - 44 - 45 - 46 - 47 - 48 - 49 - 50 - 51 - 52 - 53 - 54 - 55 - 56 - 57 - 58 - 59 - 60 - 61 - 62 - 63 - 64 - 65 - 66 - 67 - 68 - 69 - 70 - 71 - 72 - 73 - 74 - 75 - 76 - 77 - 78 - 79 - 80 - 81 - 82 - 83 - 84 - 85 - 86 - 87 - 88 - 89 - 90 - 91 - 92 - 93 - 94 - 95 - 96 - 97 - 98 - 99 - 100 - 101 - 102 - 103 - 104 - 105 - 106 - 107 - 108 - 109 - 110 - 111 - 112 - 113 - 114 - 115 - 116 - 117 - 118 - 119 - 120 - 121 - 122 >=> 123

julia> show_dynkin_diagram(:C, 2)
1 <=< 2

julia> show_dynkin_diagram(:C, 3)
1 - 2 <=< 3

julia> show_dynkin_diagram(:C, 4)
1 - 2 - 3 <=< 4

julia> show_dynkin_diagram(:C, 5)
1 - 2 - 3 - 4 <=< 5

julia> show_dynkin_diagram(:C, 10)
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 <=< 10

julia> show_dynkin_diagram(:C, 123)
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 - 10 - 11 - 12 - 13 - 14 - 15 - 16 - 17 - 18 - 19 - 20 - 21 - 22 - 23 - 24 - 25 - 26 - 27 - 28 - 29 - 30 - 31 - 32 - 33 - 34 - 35 - 36 - 37 - 38 - 39 - 40 - 41 - 42 - 43 - 44 - 45 - 46 - 47 - 48 - 49 - 50 - 51 - 52 - 53 - 54 - 55 - 56 - 57 - 58 - 59 - 60 - 61 - 62 - 63 - 64 - 65 - 66 - 67 - 68 - 69 - 70 - 71 - 72 - 73 - 74 - 75 - 76 - 77 - 78 - 79 - 80 - 81 - 82 - 83 - 84 - 85 - 86 - 87 - 88 - 89 - 90 - 91 - 92 - 93 - 94 - 95 - 96 - 97 - 98 - 99 - 100 - 101 - 102 - 103 - 104 - 105 - 106 - 107 - 108 - 109 - 110 - 111 - 112 - 113 - 114 - 115 - 116 - 117 - 118 - 119 - 120 - 121 - 122 <=< 123

julia> show_dynkin_diagram(:D, 4)
.     3
     /
1 - 2
     \
      4

julia> show_dynkin_diagram(:D, 5)
.         4
         /
1 - 2 - 3
         \
          5

julia> show_dynkin_diagram(:D, 10)
.                             9
                             /
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8
                             \
                              10

julia> show_dynkin_diagram(:D, 123)
.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       122
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       /
1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 - 10 - 11 - 12 - 13 - 14 - 15 - 16 - 17 - 18 - 19 - 20 - 21 - 22 - 23 - 24 - 25 - 26 - 27 - 28 - 29 - 30 - 31 - 32 - 33 - 34 - 35 - 36 - 37 - 38 - 39 - 40 - 41 - 42 - 43 - 44 - 45 - 46 - 47 - 48 - 49 - 50 - 51 - 52 - 53 - 54 - 55 - 56 - 57 - 58 - 59 - 60 - 61 - 62 - 63 - 64 - 65 - 66 - 67 - 68 - 69 - 70 - 71 - 72 - 73 - 74 - 75 - 76 - 77 - 78 - 79 - 80 - 81 - 82 - 83 - 84 - 85 - 86 - 87 - 88 - 89 - 90 - 91 - 92 - 93 - 94 - 95 - 96 - 97 - 98 - 99 - 100 - 101 - 102 - 103 - 104 - 105 - 106 - 107 - 108 - 109 - 110 - 111 - 112 - 113 - 114 - 115 - 116 - 117 - 118 - 119 - 120 - 121
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        123     

julia> show_dynkin_diagram(:E, 6)
1 - 3 - 4 - 5 - 6
        |
        2

julia> show_dynkin_diagram(:E, 7)
1 - 3 - 4 - 5 - 6 - 7
        |
        2

julia> show_dynkin_diagram(:E, 8)
1 - 3 - 4 - 5 - 6 - 7 - 8
        |
        2

julia> show_dynkin_diagram(:F, 4)
1 - 2 >=> 3 - 4

julia> show_dynkin_diagram(:G, 2)
1 <<< 2
```
"""
function dummy_placeholder end

end

@testset "LieAlgebras.DynkinDiagram" begin
  @testset "show and print Dynkin diagrams" begin
    # temporarily disable GC logging to avoid glitches in the doctests
    VERSION >= v"1.8.0" && GC.enable_logging(false)
    doctest(nothing, [AuxDocTest_DynkinDiagram])
    #doctest(nothing, [AuxDocTest_DynkinDiagram]; fix=true)
    VERSION >= v"1.8.0" && GC.enable_logging(true)
  end
end
