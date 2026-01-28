```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```
# Realizations of self-projecting matroids
This collection contains the (self-projecting) realization spaces of self-projecting matroids of rank k on n elements over characteristic zero for (k,n) in {(2,4),...,(2,12),(3,6),(3,7),(3,8),(4,8),(4,9),(5,10)}.
It accompanies the article "The self-projecting Grassmannian" by Alheydis Geiger and Francesca Zaffalon [GZ25](@cite).

**Warning:** The database is still under construction. The collections for selfprojecting matroids of rank 4 on 9 elements and for selfprojecting matroids of rank 5 on 10 elements are not complete yet, but underway.

For the cases {(2,4),(3,6),(4,8),(5,10)} the database stores material from the  article [GHSV24](@cite), for which the accompanying code (in Macaulay2, Magma, Matlab, OSCAR and SageMath) can be found on [github](https://github.com/sachihashimoto/self-dual).

## How to access the database
After installing OSCAR, you can access the database as described here https://docs.oscar-system.org/dev/Experimental/OscarDB/introduction/#get_db
The realization spaces of the matroids are contained in the collection
 ``Combinatorics.SelfProjectingMatroids``

You can query the database using the following parameters
 * identifier of the database entry ``data.name``
 * rank of the matroid ``data.rank``
 * size of the groundset of the matroid ``data.length_groundset``
 * dimension of its realization space ``data.dim_r``
 * dimension of its self-projecting realization space ``data.dim_s``
 * whether the realization space and the self-projecting realization space are equal ``data.equality_of_realizationspaces``
 
 Note that all query entries in the according dictionaries are strings, except if the value asked for is ``nothing``.

 Warning: for rank 3 on 8 elements, for rank 4 on 9 elements and for rank 5 on 10 elements the computation of the selfprojecting realization space did not always terminate. In these cases (as in the example above) the proeprties that could not be computed, like``dim_s``, ``equality_of_realizationspaces`` and 
``selfprojecting_realization_space``, are set to ``nothing``.

```julia-repl
julia> r3n8 = find_one(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"3", "data.length_groundset"=>"8", "data.dim_s"=>nothing]))
The matroid is of rank 3 on 8 elements.
The realization space is
  [1   0   0   1       1      1    x[4]    x[5]]
  [0   1   0   1   x[12]   x[8]       1       1]
  [0   0   1   1   x[12]      1   x[15]   x[15]]
in the multivariate polynomial ring in 15 variables over QQ
within the vanishing set of the ideal
Ideal with 10 generators
avoiding the zero loci of the polynomials
RingElem[x[12], x[15], -x[8], -x[8] + 1, x[15] - 1, x[4], x[5], -x[12] + 1, x[4] - x[15], x[5] - x[15], -x[4] + 1, -x[5] + 1, x[8] - x[12], x[4] - x[5], x[8]*x[15] - 1, x[4]*x[12] - x[15], x[5]*x[12] - x[15], -x[4]*x[12] + 1, -x[5]*x[12] + 1, -x[4]*x[8] + 1, -x[5]*x[8] + 1, -x[4]*x[8]*x[12] + x[4]*x[12] + x[8]*x[15] - x[12]*x[15] + x[12] - 1, -x[5]*x[8]*x[12] + x[5]*x[12] + x[8]*x[15] - x[12]*x[15] + x[12] - 1]
The computation of the self-projecting realization space did not terminate.
```


Once you have decided on the database entry you want to investigate more closely you have the following options.
```julia-repl
julia> MR = find_one(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"4","data.length_groundset"=>"9","data.equality_of_realizationspaces"=>"false"]))
The matroid is of rank 4 on 9 elements.
The realization space is
  [1   0   0   0   1   x[3]    x[3]       1       1]
  [0   1   0   0   1   x[9]       1    x[9]       1]
  [0   0   1   0   1      1   x[13]   x[14]   x[15]]
  [0   0   0   1   1      1       1       1       0]
in the multivariate polynomial ring in 20 variables over QQ
within the vanishing set of the ideal
Ideal with 15 generators
avoiding the zero loci of the polynomials
RingElem[-x[13], -x[14], -x[15], -x[13] + 1, -x[14] + 1, x[13] - x[14], x[9], x[9] - 1, -x[9] + x[14], x[15] - 1, -x[13] + x[15], -x[3], -x[3] + 1, x[3] - x[13], x[14] - x[15], -x[3] + x[9], -x[3] - x[9] + 2, x[9]*x[13] - 1, x[9]*x[15] - 1, -x[9]*x[13] + x[14], x[9]*x[15] - x[14], x[9]*x[15] - x[14] - x[15] + 1, -x[9]*x[15] - x[13] + x[15] + 1, x[9]*x[15] + x[13] - x[14] - x[15], -x[3]*x[14] + 1, -x[3]*x[15] + 1, -x[3]*x[14] + x[13], -x[3]*x[15] + x[13], -x[3]*x[15] + x[13] + x[15] - 1, x[3]*x[15] + x[14] - x[15] - 1, x[3]*x[15] - x[13] + x[14] - x[15], x[3]*x[9] - 1, x[3]*x[9] + x[3]*x[13] - 2*x[3] - x[9]*x[13] + 1, -x[3]*x[9] + x[3]*x[14] - x[9]*x[14] + 2*x[9] - 1, -x[3]*x[9] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + 1, x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - x[3]*x[9] - x[3]*x[14] - x[9]*x[13] + 1, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] + 1, -x[3]*x[9]*x[15] + x[3]*x[14] - x[9]*x[14] + x[9]*x[15] + x[9] - 1, -x[3]*x[9]*x[15] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + x[15], x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - 2*x[3]*x[9] - x[3]*x[13] - x[3]*x[14] + 2*x[3] - x[9]*x[13] - x[9]*x[14] + 2*x[9] + x[13] + x[14] - 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] - x[9]*x[15] + x[9] + x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] - x[9]*x[14] + x[9]*x[15] + x[9] - x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] + x[9]*x[13] + x[9]*x[15] - x[9] - x[13] - x[14] - x[15] + 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] + x[9]*x[14] - x[9]*x[15] - x[9] - x[13] - x[14] + x[15] + 2]
The selfprojecting realization space is
  [1   0   0   0   1   x[3]    x[3]       1       1]
  [0   1   0   0   1   x[9]       1    x[9]       1]
  [0   0   1   0   1      1   x[13]   x[14]   x[15]]
  [0   0   0   1   1      1       1       1       0]
in the multivariate polynomial ring in 20 variables over QQ
within the vanishing set of the ideal
Ideal with 24 generators
avoiding the zero loci of the polynomials
RingElem[-x[13], -x[14], -x[15], -x[13] + 1, -x[14] + 1, x[13] - x[14], x[9], x[9] - 1, -x[9] + x[14], x[15] - 1, -x[13] + x[15], -x[3], -x[3] + 1, x[3] - x[13], x[14] - x[15], -x[3] + x[9], -x[3] - x[9] + 2, x[9]*x[13] - 1, x[9]*x[15] - 1, -x[9]*x[13] + x[14], x[9]*x[15] - x[14], x[9]*x[15] - x[14] - x[15] + 1, -x[9]*x[15] - x[13] + x[15] + 1, x[9]*x[15] + x[13] - x[14] - x[15], -x[3]*x[14] + 1, -x[3]*x[15] + 1, -x[3]*x[14] + x[13], -x[3]*x[15] + x[13], -x[3]*x[15] + x[13] + x[15] - 1, x[3]*x[15] + x[14] - x[15] - 1, x[3]*x[15] - x[13] + x[14] - x[15], x[3]*x[9] - 1, x[3]*x[9] + x[3]*x[13] - 2*x[3] - x[9]*x[13] + 1, -x[3]*x[9] + x[3]*x[14] - x[9]*x[14] + 2*x[9] - 1, -x[3]*x[9] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + 1, x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - x[3]*x[9] - x[3]*x[14] - x[9]*x[13] + 1, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] + 1, -x[3]*x[9]*x[15] + x[3]*x[14] - x[9]*x[14] + x[9]*x[15] + x[9] - 1, -x[3]*x[9]*x[15] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + x[15], x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - 2*x[3]*x[9] - x[3]*x[13] - x[3]*x[14] + 2*x[3] - x[9]*x[13] - x[9]*x[14] + 2*x[9] + x[13] + x[14] - 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] - x[9]*x[15] + x[9] + x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] - x[9]*x[14] + x[9]*x[15] + x[9] - x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] + x[9]*x[13] + x[9]*x[15] - x[9] - x[13] - x[14] - x[15] + 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] + x[9]*x[14] - x[9]*x[15] - x[9] - x[13] - x[14] + x[15] + 2]
The closures of the realization space and the self-projecting realization space are not equal.

julia> name(MR)
"r_4_n_9_0002"

julia> Oscar.matroid(MR)
Matroid of rank 4 on 9 elements

julia> rank(MR)
4

julia> length_groundset(MR)
9

julia> realization_space(MR)
The realization space is
  [1   0   0   0   1   x[3]    x[3]       1       1]
  [0   1   0   0   1   x[9]       1    x[9]       1]
  [0   0   1   0   1      1   x[13]   x[14]   x[15]]
  [0   0   0   1   1      1       1       1       0]
in the multivariate polynomial ring in 20 variables over QQ
within the vanishing set of the ideal
Ideal with 15 generators
avoiding the zero loci of the polynomials
RingElem[-x[13], -x[14], -x[15], -x[13] + 1, -x[14] + 1, x[13] - x[14], x[9], x[9] - 1, -x[9] + x[14], x[15] - 1, -x[13] + x[15], -x[3], -x[3] + 1, x[3] - x[13], x[14] - x[15], -x[3] + x[9], -x[3] - x[9] + 2, x[9]*x[13] - 1, x[9]*x[15] - 1, -x[9]*x[13] + x[14], x[9]*x[15] - x[14], x[9]*x[15] - x[14] - x[15] + 1, -x[9]*x[15] - x[13] + x[15] + 1, x[9]*x[15] + x[13] - x[14] - x[15], -x[3]*x[14] + 1, -x[3]*x[15] + 1, -x[3]*x[14] + x[13], -x[3]*x[15] + x[13], -x[3]*x[15] + x[13] + x[15] - 1, x[3]*x[15] + x[14] - x[15] - 1, x[3]*x[15] - x[13] + x[14] - x[15], x[3]*x[9] - 1, x[3]*x[9] + x[3]*x[13] - 2*x[3] - x[9]*x[13] + 1, -x[3]*x[9] + x[3]*x[14] - x[9]*x[14] + 2*x[9] - 1, -x[3]*x[9] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + 1, x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - x[3]*x[9] - x[3]*x[14] - x[9]*x[13] + 1, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] + 1, -x[3]*x[9]*x[15] + x[3]*x[14] - x[9]*x[14] + x[9]*x[15] + x[9] - 1, -x[3]*x[9]*x[15] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + x[15], x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - 2*x[3]*x[9] - x[3]*x[13] - x[3]*x[14] + 2*x[3] - x[9]*x[13] - x[9]*x[14] + 2*x[9] + x[13] + x[14] - 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] - x[9]*x[15] + x[9] + x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] - x[9]*x[14] + x[9]*x[15] + x[9] - x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] + x[9]*x[13] + x[9]*x[15] - x[9] - x[13] - x[14] - x[15] + 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] + x[9]*x[14] - x[9]*x[15] - x[9] - x[13] - x[14] + x[15] + 2]

julia> dim_r(MR)
5

julia> selfprojecting_realization_space(MR)
The selfprojecting realization space is
  [1   0   0   0   1   x[3]    x[3]       1       1]
  [0   1   0   0   1   x[9]       1    x[9]       1]
  [0   0   1   0   1      1   x[13]   x[14]   x[15]]
  [0   0   0   1   1      1       1       1       0]
in the multivariate polynomial ring in 20 variables over QQ
within the vanishing set of the ideal
Ideal with 24 generators
avoiding the zero loci of the polynomials
RingElem[-x[13], -x[14], -x[15], -x[13] + 1, -x[14] + 1, x[13] - x[14], x[9], x[9] - 1, -x[9] + x[14], x[15] - 1, -x[13] + x[15], -x[3], -x[3] + 1, x[3] - x[13], x[14] - x[15], -x[3] + x[9], -x[3] - x[9] + 2, x[9]*x[13] - 1, x[9]*x[15] - 1, -x[9]*x[13] + x[14], x[9]*x[15] - x[14], x[9]*x[15] - x[14] - x[15] + 1, -x[9]*x[15] - x[13] + x[15] + 1, x[9]*x[15] + x[13] - x[14] - x[15], -x[3]*x[14] + 1, -x[3]*x[15] + 1, -x[3]*x[14] + x[13], -x[3]*x[15] + x[13], -x[3]*x[15] + x[13] + x[15] - 1, x[3]*x[15] + x[14] - x[15] - 1, x[3]*x[15] - x[13] + x[14] - x[15], x[3]*x[9] - 1, x[3]*x[9] + x[3]*x[13] - 2*x[3] - x[9]*x[13] + 1, -x[3]*x[9] + x[3]*x[14] - x[9]*x[14] + 2*x[9] - 1, -x[3]*x[9] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + 1, x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - x[3]*x[9] - x[3]*x[14] - x[9]*x[13] + 1, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] + 1, -x[3]*x[9]*x[15] + x[3]*x[14] - x[9]*x[14] + x[9]*x[15] + x[9] - 1, -x[3]*x[9]*x[15] + x[3]*x[14] + x[9]*x[13] - x[13] - x[14] + x[15], x[3]*x[9]*x[13] + x[3]*x[9]*x[14] - 2*x[3]*x[9] - x[3]*x[13] - x[3]*x[14] + 2*x[3] - x[9]*x[13] - x[9]*x[14] + 2*x[9] + x[13] + x[14] - 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] - x[9]*x[13] - x[9]*x[15] + x[9] + x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] - x[9]*x[14] + x[9]*x[15] + x[9] - x[15], -x[3]*x[9]*x[15] + x[3]*x[14] + x[3]*x[15] - x[3] + x[9]*x[13] + x[9]*x[15] - x[9] - x[13] - x[14] - x[15] + 2, x[3]*x[9]*x[15] + x[3]*x[13] - x[3]*x[15] - x[3] + x[9]*x[14] - x[9]*x[15] - x[9] - x[13] - x[14] + x[15] + 2]

julia> dim_s(MR)
3

julia> equality_of_realizationspaces(MR)
false
```
The realization space $\mathcal{R}$ obtained by ``realization_space(MR)`` and the self-projecting realization space $\mathcal{S}$ obtained by ``selfprojecting_realization_space(MR)`` can be investigated using the code in the experimental section of OSCAR on MatroidRealizationSpaces.

```julia-repl
julia> R = realization_space(MR);

julia> defining_ideal(R)
Ideal generated by
  x[1] - 1
  x[2] - x[3]
  x[4] - 1
  x[5] - 1
  x[6] - 1
  x[7] - x[9]
  x[8] - 1
  x[10] - 1
  x[11] - 1
  x[12] - 1
  x[16] - 1
  x[17] - 1
  x[18] - 1
  x[19] - 1
  x[20]

julia> ambient_ring(R)
Multivariate polynomial ring in 20 variables x[1], x[2], x[3], x[4], ..., x[20]
  over rational field
  
julia> S = selfprojecting_realization_space(MR);

julia> defining_ideal(S)
Ideal generated by
  x[1] - 1
  x[2] - x[3]
  x[3]^2*x[9]*x[13] - x[3]^2*x[9] - x[3]^2*x[13]*x[14] + x[3]^2*x[14] - x[3]*x[9
  ]^2*x[14] + x[3]*x[9]^2 - 2*x[3]*x[9]*x[13] + 2*x[3]*x[9]*x[14] + 2*x[3]*x[13]
  *x[14] - 3*x[3]*x[14] + x[3] + x[9]^2*x[13]*x[14] - x[9]^2*x[13] - 2*x[9]*x[13
  ]*x[14] + 3*x[9]*x[13] - x[9] - x[13] + x[14]
  x[3]^2*x[9]*x[15] - x[3]^2*x[14]*x[15] + x[3]*x[9]*x[14] - x[3]*x[9]*x[15]^2 -
   x[3]*x[9]*x[15] - x[3]*x[9] - x[3]*x[13]*x[14]*x[15] + x[3]*x[13]*x[14] + x[3
  ]*x[14]*x[15]^2 + x[3]*x[14]*x[15] - x[3]*x[14] + x[3]*x[15] - x[9]*x[13]*x[14
  ]*x[15] + x[9]*x[13]*x[15]^2 + x[9]*x[15] + x[13]*x[14] - x[13] - x[14] - x[15
  ]^2 + 1
  x[3]^2*x[13]*x[14]*x[15] - x[3]^2*x[13]*x[14] - x[3]^2*x[14]^2*x[15] + x[3]^2*
  x[14]*x[15] + x[3]^2*x[14] - x[3]^2*x[15] + x[3]*x[9]*x[14]^2 - x[3]*x[9]*x[14
  ] - x[3]*x[9]*x[15]^2 + x[3]*x[9]*x[15] + x[3]*x[13]*x[14]^2 - 3*x[3]*x[13]*x[
  14]*x[15] + x[3]*x[13]*x[14] + x[3]*x[13] + x[3]*x[14]^2*x[15] - x[3]*x[14]^2
  + x[3]*x[14]*x[15]^2 - 2*x[3]*x[14] + x[3]*x[15] - x[9]*x[13]*x[14]^2 + x[9]*x
  [13]*x[14] + x[9]*x[13]*x[15]^2 - x[9]*x[13]*x[15] + x[13]*x[14] + x[13]*x[15]
   - 2*x[13] - x[14]*x[15] + x[14] - x[15]^2 + x[15]
  x[3]*x[9]^2*x[15] + x[3]*x[9]*x[13] - x[3]*x[9]*x[15]^2 - x[3]*x[9]*x[15] - x[
  3]*x[9] - x[3]*x[13]*x[14]*x[15] + x[3]*x[14]*x[15]^2 + x[3]*x[15] - x[9]^2*x[
  13]*x[15] - x[9]*x[13]*x[14]*x[15] + x[9]*x[13]*x[14] + x[9]*x[13]*x[15]^2 + x
  [9]*x[13]*x[15] - x[9]*x[13] + x[9]*x[15] + x[13]*x[14] - x[13] - x[14] - x[15
  ]^2 + 1
  x[3]*x[9]*x[13]^2 - x[3]*x[9]*x[13] - x[3]*x[9]*x[15]^2 + x[3]*x[9]*x[15] - x[
  3]*x[13]^2*x[14] + x[3]*x[13]*x[14] + x[3]*x[14]*x[15]^2 - x[3]*x[14]*x[15] -
  x[9]^2*x[13]^2*x[15] + x[9]^2*x[13]*x[14]*x[15] - x[9]^2*x[13]*x[14] + x[9]^2*
  x[13]*x[15] + x[9]^2*x[13] - x[9]^2*x[15] + x[9]*x[13]^2*x[14] + x[9]*x[13]^2*
  x[15] - x[9]*x[13]^2 - 3*x[9]*x[13]*x[14]*x[15] + x[9]*x[13]*x[14] + x[9]*x[13
  ]*x[15]^2 - 2*x[9]*x[13] + x[9]*x[14] + x[9]*x[15] + x[13]*x[14] - x[13]*x[15]
   + x[13] + x[14]*x[15] - 2*x[14] - x[15]^2 + x[15]
  x[3]*x[9]*x[13]*x[15] - x[3]*x[9]*x[15]^2 - x[3]*x[13]*x[14]*x[15] + x[3]*x[14
  ]*x[15]^2 - x[9]*x[13]*x[14]*x[15] + x[9]*x[13]*x[14] + x[9]*x[13]*x[15]^2 - x
  [9]*x[13]*x[15] - x[9]*x[13] + x[9]*x[15] + x[13]*x[14] - x[14] - x[15]^2 + x[
  15]
  x[3]*x[9]*x[14]*x[15] - x[3]*x[9]*x[15]^2 - x[3]*x[13]*x[14]*x[15] + x[3]*x[13
  ]*x[14] + x[3]*x[14]*x[15]^2 - x[3]*x[14]*x[15] - x[3]*x[14] + x[3]*x[15] - x[
  9]*x[13]*x[14]*x[15] + x[9]*x[13]*x[15]^2 + x[13]*x[14] - x[13] - x[15]^2 + x[
  15]
  x[3]*x[13]^2*x[14]*x[15] - x[3]*x[13]^2*x[14] - x[3]*x[13]*x[14]^2*x[15] - x[3
  ]*x[13]*x[14]*x[15]^2 + 2*x[3]*x[13]*x[14]*x[15] + x[3]*x[13]*x[14] - x[3]*x[1
  3]*x[15] + x[3]*x[14]^2*x[15]^2 - x[3]*x[14]*x[15]^2 - x[3]*x[14]*x[15] + x[3]
  *x[15]^2 + x[9]*x[13]^2*x[14]*x[15] - x[9]*x[13]^2*x[15]^2 - x[9]*x[13]*x[14]^
  2*x[15] + x[9]*x[13]*x[14]^2 + x[9]*x[13]*x[14]*x[15]^2 - 2*x[9]*x[13]*x[14]*x
  [15] - x[9]*x[13]*x[14] + x[9]*x[13]*x[15]^2 + x[9]*x[13]*x[15] + x[9]*x[14]*x
  [15] - x[9]*x[15]^2 - x[13]^2*x[14] + x[13]^2 + x[13]*x[14]^2 + x[13]*x[15]^2
  - 2*x[13]*x[15] - x[14]^2 - x[14]*x[15]^2 + 2*x[14]*x[15]
  x[4] - 1
  x[5] - 1
  x[6] - 1
  x[7] - x[9]
  x[8] - 1
  x[9]^2*x[13]^2*x[15]^2 - x[9]^2*x[13]*x[14]*x[15]^2 + x[9]^2*x[13]*x[14]*x[15]
   - x[9]^2*x[13]*x[15]^2 - x[9]^2*x[13]*x[15] + x[9]^2*x[15]^2 - 2*x[9]*x[13]^2
  *x[14]*x[15] + x[9]*x[13]^2*x[14] - x[9]*x[13]^2 + 2*x[9]*x[13]*x[14]*x[15]^2
  + x[9]*x[13]*x[14]*x[15] - x[9]*x[13]*x[14] - 2*x[9]*x[13]*x[15]^2 + 3*x[9]*x[
  13]*x[15] + x[9]*x[13] - x[9]*x[14]*x[15] - x[9]*x[15] + x[13]^2*x[14] - 2*x[1
  3]*x[14] - x[14]*x[15]^2 + x[14]*x[15] + x[14] + x[15]^2 - x[15]
  x[10] - 1
  x[11] - 1
  x[12] - 1
  x[16] - 1
  x[17] - 1
  x[18] - 1
  x[19] - 1
  x[20]
```
The example above showcased $\mathcal{R}\supsetneq\mathcal{S}$. Below you see an example with $\mathcal{R}=\mathcal{S}$.
```julia-repl
julia> MR = find_one(db["Combinatorics.SelfProjectingMatroids"], Dict("data.name"=>"r_3_n_8_10"))
The matroid is of rank 3 on 8 elements.
The realization space is
  [1   0   0   1       1   x[14]    x[4]       1]
  [0   1   0   1   x[12]       1       1       1]
  [0   0   1   1   x[12]   x[14]   x[14]   x[15]]
in the multivariate polynomial ring in 15 variables over QQ
within the vanishing set of the ideal
Ideal with 11 generators
avoiding the zero loci of the polynomials
RingElem[x[12], x[14], x[15], x[14] - 1, x[15] - 1, -x[14] + x[15], x[4], -x[12] + 1, x[4] - x[14], x[12] - x[15], -x[4] + 1, x[4]*x[12] - x[14], -x[4]*x[15] + x[14], -x[12]*x[14] + 1, -x[4]*x[12] + 1, -x[12]*x[14]*x[15] + 2*x[12]*x[14] - x[12] - x[14] + x[15], -x[4]*x[12]*x[15] + x[4]*x[12] + x[12]*x[14] - x[12] - x[14] + x[15]]
The selfprojecting realization space is
  [1   0   0   1       1   x[14]    x[4]       1]
  [0   1   0   1   x[12]       1       1       1]
  [0   0   1   1   x[12]   x[14]   x[14]   x[15]]
in the multivariate polynomial ring in 15 variables over QQ
within the vanishing set of the ideal
Ideal with 11 generators
avoiding the zero loci of the polynomials
RingElem[x[12], x[14], x[15], x[14] - 1, x[15] - 1, -x[14] + x[15], x[4], -x[12] + 1, x[4] - x[14], x[12] - x[15], -x[4] + 1, x[4]*x[12] - x[14], -x[4]*x[15] + x[14], -x[12]*x[14] + 1, -x[4]*x[12] + 1, -x[12]*x[14]*x[15] + 2*x[12]*x[14] - x[12] - x[14] + x[15], -x[4]*x[12]*x[15] + x[4]*x[12] + x[12]*x[14] - x[12] - x[14] + x[15]]
The closures of the realization space and the self-projecting realization space are equal.

julia> name(MR)
"r_3_n_8_10"

julia> Oscar.matroid(MR)
Matroid of rank 3 on 8 elements

julia> rank(MR)
3

julia> length_groundset(MR)
8

julia> realization_space(MR)
The realization space is
  [1   0   0   1       1   x[14]    x[4]       1]
  [0   1   0   1   x[12]       1       1       1]
  [0   0   1   1   x[12]   x[14]   x[14]   x[15]]
in the multivariate polynomial ring in 15 variables over QQ
within the vanishing set of the ideal
Ideal with 11 generators
avoiding the zero loci of the polynomials
RingElem[x[12], x[14], x[15], x[14] - 1, x[15] - 1, -x[14] + x[15], x[4], -x[12] + 1, x[4] - x[14], x[12] - x[15], -x[4] + 1, x[4]*x[12] - x[14], -x[4]*x[15] + x[14], -x[12]*x[14] + 1, -x[4]*x[12] + 1, -x[12]*x[14]*x[15] + 2*x[12]*x[14] - x[12] - x[14] + x[15], -x[4]*x[12]*x[15] + x[4]*x[12] + x[12]*x[14] - x[12] - x[14] + x[15]]

julia> dim_r(MR)
4

julia> selfprojecting_realization_space(MR)
The selfprojecting realization space is
  [1   0   0   1       1   x[14]    x[4]       1]
  [0   1   0   1   x[12]       1       1       1]
  [0   0   1   1   x[12]   x[14]   x[14]   x[15]]
in the multivariate polynomial ring in 15 variables over QQ
within the vanishing set of the ideal
Ideal with 11 generators
avoiding the zero loci of the polynomials
RingElem[x[12], x[14], x[15], x[14] - 1, x[15] - 1, -x[14] + x[15], x[4], -x[12] + 1, x[4] - x[14], x[12] - x[15], -x[4] + 1, x[4]*x[12] - x[14], -x[4]*x[15] + x[14], -x[12]*x[14] + 1, -x[4]*x[12] + 1, -x[12]*x[14]*x[15] + 2*x[12]*x[14] - x[12] - x[14] + x[15], -x[4]*x[12]*x[15] + x[4]*x[12] + x[12]*x[14] - x[12] - x[14] + x[15]]

julia> dim_s(MR)
4

julia> equality_of_realizationspaces(MR)
true
```

## How to verify claims from the article
To verify Tables 2, 3, and 4 from the article, one can use queries to the database. The example below shows how to generate the line of Table 2 with respect to the matroids of rank 3 on 7 elements. Recall that the uniform matroids are not stored in the database.
The other rows as well as Table 3 can be verified similarly. Note that the database collection for (4,9) is not filled completely yet.
```julia-repl
julia> t2 = [length(db["Combinatorics.SelfProjectingMatroids"], Dict("data.rank"=>"3", "data.length_groundset"=>"7","data.dim_s"=>"$i")) for i in -1:5]
7-element Vector{Int64}:
 1
 1
 1
 3
 3
 1
 1
```
To obtain numbers from Table 4, i.e. the distribution of realizable matroids without selfprojecting realization you can use the following code:
```julia-repl
julia> t4 = [length(db["Combinatorics.SelfProjectingMatroids"], Dict("data.rank"=>"4", "data.length_groundset"=>"9","data.dim_r"=>"$i","data.dim_s"=>"-1")) for i in 0:6]
7-element Vector{Int64}:
    4
  103
  494
 1089
  738
  124
    4
```
One can count the matroids for which $\mathcal{R}$ and $\mathcal{S}$ are known and do not coincide as follows:
```julia-repl
julia> length(db["Combinatorics.SelfProjectingMatroids"], Dict("data.rank"=>"4", "data.length_groundset"=>"9","data.equality_of_realizationspaces"=>"false"))
5399
```
Theorem 4.11 claims that there are at least 5400 matroids of rank 4 on 9 elements with $\mathcal{R}\supsetneq\mathcal{S}$. Since the database does not count the uniform matroid $U_{4,9}$, the claim is verified.

The user can find the selfprojecting matroids for which the computation of the selfprojecting realization space was too costly and did not terminate as follows:
```julia-repl
julia> notterminated = length([r for r in find(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"3", "data.length_groundset"=>"8","data.dim_s"=>nothing]))])
4
```

**Example 4.12**

In order to work with the database and/or compute self-projecting realization spaces of matroids in OSCAR, you need to use version 1.6.0 or later.
To reproduce example 4.12 you can access the relevant file from the database.
```julia-repl
julia> using Oscar

julia> db = Oscar.OscarDB.get_db();

julia> find_one(db["Combinatorics.SelfProjectingMatroids"], Dict(["name"=>"r_4_n_9_index_5985"]))
The matroid is of rank 4 on 9 elements.
The realization space is
  [1   0   0   0   2//3   0      1   1   1//2]
  [0   1   0   0      0   2   1//2   1   1//2]
  [0   0   1   0      1   1      1   1      1]
  [0   0   0   1      2   2      2   1      1]
in the multivariate polynomial ring in 20 variables over QQ
within the vanishing set of the ideal
Ideal with 20 generators
avoiding the zero loci of the polynomials
RingElem[2]
The matroid does not have a self-projecting realization over characteristic zero.
The closures of the realization space and the self-projecting realization space are not equal.
```

## Additional Code

Additional code as well as the original code and output from the computations in magma can be found at https://github.com/AlheydisGeiger/selfprojectingGrassmannian

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Alheydis Geiger](https://www.mis.mpg.de/people/alheydis-geiger),
* [Francesca Zaffalon](https://sites.google.com/view/francescazaffalon/home).


You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).


Software used to create the database collection ``Combinatorics.SelfProjectingMatroids``: Magma (V2.27), Julia (Version 1.12.1), OSCAR (version 1.6.0-DEV), 
GNU parallel 20221122


Last updated 11/12/2025.

