@testset "LieAlgebras.DynkinDiagram" begin

#! format: off
Oscar.@_AuxDocTest "show and print Dynkin diagrams", (fix = false),
raw"""
canonical labels

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

non-canonical labels

```jldoctest show_dynkin_diagram_with_labels.test
julia> using Oscar

julia> show_dynkin_diagram(:A, 1, [5])
5

julia> show_dynkin_diagram(:A, 2, [42, 17])
42 - 17

julia> show_dynkin_diagram(:A, 3, [3,2,1])
3 - 2 - 1

julia> show_dynkin_diagram(:A, 4, [3,12,101,2])
3 - 12 - 101 - 2

julia> show_dynkin_diagram(:A, 5, 100:44:300)
100 - 144 - 188 - 232 - 276

julia> show_dynkin_diagram(:A, 10, 100:-10:10)
100 - 90 - 80 - 70 - 60 - 50 - 40 - 30 - 20 - 10

julia> show_dynkin_diagram(:A, 123, 122:-2:-122)
122 - 120 - 118 - 116 - 114 - 112 - 110 - 108 - 106 - 104 - 102 - 100 - 98 - 96 - 94 - 92 - 90 - 88 - 86 - 84 - 82 - 80 - 78 - 76 - 74 - 72 - 70 - 68 - 66 - 64 - 62 - 60 - 58 - 56 - 54 - 52 - 50 - 48 - 46 - 44 - 42 - 40 - 38 - 36 - 34 - 32 - 30 - 28 - 26 - 24 - 22 - 20 - 18 - 16 - 14 - 12 - 10 - 8 - 6 - 4 - 2 - 0 - -2 - -4 - -6 - -8 - -10 - -12 - -14 - -16 - -18 - -20 - -22 - -24 - -26 - -28 - -30 - -32 - -34 - -36 - -38 - -40 - -42 - -44 - -46 - -48 - -50 - -52 - -54 - -56 - -58 - -60 - -62 - -64 - -66 - -68 - -70 - -72 - -74 - -76 - -78 - -80 - -82 - -84 - -86 - -88 - -90 - -92 - -94 - -96 - -98 - -100 - -102 - -104 - -106 - -108 - -110 - -112 - -114 - -116 - -118 - -120 - -122

julia> show_dynkin_diagram(:B, 2, [42, 17])
42 >=> 17

julia> show_dynkin_diagram(:B, 3, [3,2,1])
3 - 2 >=> 1

julia> show_dynkin_diagram(:B, 4, [3,12,101,2])
3 - 12 - 101 >=> 2

julia> show_dynkin_diagram(:B, 5, 100:44:300)
100 - 144 - 188 - 232 >=> 276

julia> show_dynkin_diagram(:B, 10, 100:-10:10)
100 - 90 - 80 - 70 - 60 - 50 - 40 - 30 - 20 >=> 10

julia> show_dynkin_diagram(:B, 123, 122:-2:-122)
122 - 120 - 118 - 116 - 114 - 112 - 110 - 108 - 106 - 104 - 102 - 100 - 98 - 96 - 94 - 92 - 90 - 88 - 86 - 84 - 82 - 80 - 78 - 76 - 74 - 72 - 70 - 68 - 66 - 64 - 62 - 60 - 58 - 56 - 54 - 52 - 50 - 48 - 46 - 44 - 42 - 40 - 38 - 36 - 34 - 32 - 30 - 28 - 26 - 24 - 22 - 20 - 18 - 16 - 14 - 12 - 10 - 8 - 6 - 4 - 2 - 0 - -2 - -4 - -6 - -8 - -10 - -12 - -14 - -16 - -18 - -20 - -22 - -24 - -26 - -28 - -30 - -32 - -34 - -36 - -38 - -40 - -42 - -44 - -46 - -48 - -50 - -52 - -54 - -56 - -58 - -60 - -62 - -64 - -66 - -68 - -70 - -72 - -74 - -76 - -78 - -80 - -82 - -84 - -86 - -88 - -90 - -92 - -94 - -96 - -98 - -100 - -102 - -104 - -106 - -108 - -110 - -112 - -114 - -116 - -118 - -120 >=> -122

julia> show_dynkin_diagram(:C, 2, [42, 17])
42 <=< 17

julia> show_dynkin_diagram(:C, 3, [3,2,1])
3 - 2 <=< 1

julia> show_dynkin_diagram(:C, 4, [3,12,101,2])
3 - 12 - 101 <=< 2

julia> show_dynkin_diagram(:C, 5, 100:44:300)
100 - 144 - 188 - 232 <=< 276

julia> show_dynkin_diagram(:C, 10, 100:-10:10)
100 - 90 - 80 - 70 - 60 - 50 - 40 - 30 - 20 <=< 10

julia> show_dynkin_diagram(:C, 123, 122:-2:-122)
122 - 120 - 118 - 116 - 114 - 112 - 110 - 108 - 106 - 104 - 102 - 100 - 98 - 96 - 94 - 92 - 90 - 88 - 86 - 84 - 82 - 80 - 78 - 76 - 74 - 72 - 70 - 68 - 66 - 64 - 62 - 60 - 58 - 56 - 54 - 52 - 50 - 48 - 46 - 44 - 42 - 40 - 38 - 36 - 34 - 32 - 30 - 28 - 26 - 24 - 22 - 20 - 18 - 16 - 14 - 12 - 10 - 8 - 6 - 4 - 2 - 0 - -2 - -4 - -6 - -8 - -10 - -12 - -14 - -16 - -18 - -20 - -22 - -24 - -26 - -28 - -30 - -32 - -34 - -36 - -38 - -40 - -42 - -44 - -46 - -48 - -50 - -52 - -54 - -56 - -58 - -60 - -62 - -64 - -66 - -68 - -70 - -72 - -74 - -76 - -78 - -80 - -82 - -84 - -86 - -88 - -90 - -92 - -94 - -96 - -98 - -100 - -102 - -104 - -106 - -108 - -110 - -112 - -114 - -116 - -118 - -120 <=< -122

julia> show_dynkin_diagram(:D, 4, [3,12,101,2])
.      101
      /
3 - 12
      \
       2

julia> show_dynkin_diagram(:D, 5, 100:44:300)
.               232
               /
100 - 144 - 188
               \
                276

julia> show_dynkin_diagram(:D, 10, 100:-10:10)
.                                      20
                                      /
100 - 90 - 80 - 70 - 60 - 50 - 40 - 30
                                      \
                                       10

julia> show_dynkin_diagram(:D, 123, 122:-2:-122)
.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  -120
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  /
122 - 120 - 118 - 116 - 114 - 112 - 110 - 108 - 106 - 104 - 102 - 100 - 98 - 96 - 94 - 92 - 90 - 88 - 86 - 84 - 82 - 80 - 78 - 76 - 74 - 72 - 70 - 68 - 66 - 64 - 62 - 60 - 58 - 56 - 54 - 52 - 50 - 48 - 46 - 44 - 42 - 40 - 38 - 36 - 34 - 32 - 30 - 28 - 26 - 24 - 22 - 20 - 18 - 16 - 14 - 12 - 10 - 8 - 6 - 4 - 2 - 0 - -2 - -4 - -6 - -8 - -10 - -12 - -14 - -16 - -18 - -20 - -22 - -24 - -26 - -28 - -30 - -32 - -34 - -36 - -38 - -40 - -42 - -44 - -46 - -48 - -50 - -52 - -54 - -56 - -58 - -60 - -62 - -64 - -66 - -68 - -70 - -72 - -74 - -76 - -78 - -80 - -82 - -84 - -86 - -88 - -90 - -92 - -94 - -96 - -98 - -100 - -102 - -104 - -106 - -108 - -110 - -112 - -114 - -116 - -118
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   -122

julia> show_dynkin_diagram(:E, 6, [12,4,3,101,0,4])
12 - 3 - 101 - 0 - 4
         |
         4


julia> show_dynkin_diagram(:E, 7, [88,12,4,3,0,101,4])
88 - 4 - 3 - 0 - 101 - 4
         |
         12


julia> show_dynkin_diagram(:E, 8, [9,101,12,4,3,0,101,4])
9 - 12 - 4 - 3 - 0 - 101 - 4
         |
         101


julia> show_dynkin_diagram(:F, 4, [105, 2, 99, 300])
105 - 2 >=> 99 - 300

julia> show_dynkin_diagram(:G, 2, [412, 5])
412 <<< 5
```

non-simple diagram with canonical labels

```jldoctest show_ns_dynkin_diagram.test
julia> using Oscar

julia> show_dynkin_diagram([(:B, 5), (:B, 5)])
1 - 2 - 3 - 4 >=> 5

6 - 7 - 8 - 9 >=> 10

julia> show_dynkin_diagram([(:A, 2), (:B, 3), (:C, 4), (:D, 5), (:E, 6), (:F, 4), (:G, 2)])
1 - 2

3 - 4 >=> 5

6 - 7 - 8 <=< 9

.            13
            /
10 - 11 - 12
            \
             14

15 - 17 - 18 - 19 - 20
          |
          16

21 - 22 >=> 23 - 24

25 <<< 26
```

non-simple diagram with non-canonical labels

```jldoctest show_ns_dynkin_diagram_with_labels.test
julia> using Oscar

julia> show_dynkin_diagram([(:B, 5), (:B, 5)], [15,3,6,0,1000,23,8,22,65,1])
15 - 3 - 6 - 0 >=> 1000

23 - 8 - 22 - 65 >=> 1

julia> show_dynkin_diagram([(:A, 2), (:B, 3), (:C, 4), (:D, 5), (:E, 6), (:F, 4), (:G, 2)], 2*(2+3+4+5+6+4+2):-2:2)
52 - 50

48 - 46 >=> 44

42 - 40 - 38 <=< 36

.            28
            /
34 - 32 - 30
            \
             26

24 - 20 - 18 - 16 - 14
          |
          22

12 - 10 >=> 8 - 6

4 <<< 2
```

Dynkin diagram from cartan matrix

```jldoctest show_dynkin_diagram_cartan_matrix.test
julia> using Oscar

julia> show_dynkin_diagram(cartan_matrix(:A, 3))
1 - 2 - 3

julia> show_dynkin_diagram(cartan_matrix(:B, 4))
1 - 2 - 3 >=> 4

julia> show_dynkin_diagram(cartan_matrix(:C, 5))
1 - 2 - 3 - 4 <=< 5

julia> show_dynkin_diagram(cartan_matrix(:D, 6))
.             5
             /
1 - 2 - 3 - 4
             \
              6

julia> show_dynkin_diagram(cartan_matrix(:E, 7))
1 - 3 - 4 - 5 - 6 - 7
        |
        2

julia> show_dynkin_diagram(cartan_matrix(:F, 4))
1 - 2 >=> 3 - 4

julia> show_dynkin_diagram(cartan_matrix(:G, 2))
1 <<< 2

julia> show_dynkin_diagram(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2])
1 >=> 3

2 <=< 4

julia> show_dynkin_diagram(transpose(cartan_matrix(:F, 4)))
4 - 3 >=> 2 - 1

julia> show_dynkin_diagram(transpose(cartan_matrix(:G, 2)))
2 <<< 1
```

Dynkin diagram from root system

```jldoctest show_dynkin_diagram_cartan_matrix.test
julia> using Oscar

julia> show_dynkin_diagram(root_system(:A, 3))
1 - 2 - 3

julia> show_dynkin_diagram(root_system([(:B, 4)]))
1 - 2 - 3 >=> 4

julia> show_dynkin_diagram(root_system(cartan_matrix(:C, 5)))
1 - 2 - 3 - 4 <=< 5

julia> show_dynkin_diagram(root_system([(:D, 8), (:E, 6), (:F, 4), (:G, 2)]))
.                     7
                     /
1 - 2 - 3 - 4 - 5 - 6
                     \
                      8

9 - 11 - 12 - 13 - 14
         |
         10

15 - 16 >=> 17 - 18

19 <<< 20

julia> show_dynkin_diagram(root_system(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2]))
1 >=> 3

2 <=< 4

julia> show_dynkin_diagram(root_system(transpose(cartan_matrix(:F, 4))))
4 - 3 >=> 2 - 1

julia> show_dynkin_diagram(root_system(transpose(cartan_matrix(:G, 2))))
2 <<< 1
```
"""
#! format: on

end
