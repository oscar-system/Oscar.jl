```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```
# Realizations of self-projecting matroids
This page accompanies the article "The self-projecting Grassmannian" by Alheydis Geiger and Francesca Zaffalon (https://arxiv.org/abs/2511.21442). It contains the database of (self-projecting) realization spaces of self-projecting matroids of rank k on n elements over characteristic zero for (k,n) in {(2,4),...,(2,12),(3,6),(3,7),(3,8),(4,8),(4,9),(5,10)}.



**Abstract:**  We introduce the self-projecting Grassmannian, an irreducible subvariety of the Grassmannian parametrizing linear subspaces that satisfy a generalized self-duality condition. We study its relation to classical moduli spaces, such as the moduli spaces of pointed curves of genus $g$, as well as to other natural subvarieties of the Grassmannian. We further translate the self-projectivity condition into the combinatorial language of matroids, introducing self-projecting matroids, and we computationally investigate their realization spaces inside the self-projecting Grassmannian.
arXiv: <https://arxiv.org/abs/2511.21442>

**Warning:** The database is still under construction. The collections for selfprojecting matroids of rank 4 on 9 elements and for selfprojecting matroids of rank 5 on 10 elements are not complete yet.

For the cases {(2,4),(3,6),(4,8),(5,10)} the database stores material from the  article:
Alheydis Geiger, Sachi Hashimoto, Bernd Sturmfels, Raluca Vlad: Self-dual matroids from canonical curves
In: Experimental mathematics, 33 (2024) 4, p. 701-722
DOI: `10.1080/10586458.2023.2239282 <https://dx.doi.org/10.1080/10586458.2023.2239282>`_ ARXIV: https://arxiv.org/abs/2212.05910 CODE: https://github.com/sachihashimoto/self-dual
Project page created: 18/11/2025.

## How to access the database
After installing OSCAR, you can access the database as follows
```
julia> db = Oscar.OscarDB.get_db();

julia> Oscar.OscarDB.get_collection_names(db)
6-element Vector{String}:
 "SmallTreeModels"
 "Combinatorics.SelfProjectingMatroids"
 "zzlattices"
 "LeechPairs"
 "Surfaces"
 "TransitiveSimplicialComplexes"
```
You can query the database using the following parameters
 * identifier of the database entry ``data.name``
 * rank of the matroid ``data.rank``
 * size of the groundset of the matroid ``data.length_groundset``
 * dimension of its realization space ``data.dim_r``
 * dimension of its self-projecting realization space ``data.dim_s``
 * whether the realization space and the self-projecting realization space are equal ``data.equality_of_realizationspaces``
 
 Note that all query entries in the according dictionaries are strings, except if the value asked for is ``nothing``.

```
julia> r4n8 = find(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"4", "data.length_groundset"=>"8"]));
julia> length([MR for MR in r4n8])
12
```
Once you have decided on the database entry you want to investigate more closely you have the following options.
```
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
The realization space $\mathcal{R}$ obtained by ``realization_space(MR)`` and the self-projecting realization space $\mathcal{S}$ obtained by ``selfprojecting_realization_space(MR)`` can be investigated using the code in the experimental section of OSCAR on MatroidRealizationSpaces.

```
julia> R = realization_space(MR);
julia> defining_ideal(R)
Ideal generated by
  x[1] - 1
  x[2] - 1
  x[3] - x[14]
  x[5] - 1
  x[6] - 1
  x[7] - x[12]
  x[8] - 1
  x[9] - 1
  x[10] - 1
  x[11] - 1
  x[13] - x[14]
  julia> ambient_ring(R)
Multivariate polynomial ring in 15 variables x[1], x[2], x[3], x[4], ..., x[15]
  over rational field
  
julia> R = selfprojecting_realization_space(MR);
julia> defining_ideal(S)
Ideal generated by
  x[1] - 1
  x[2] - 1
  x[3] - x[14]
  x[5] - 1
  x[6] - 1
  x[7] - x[12]
  x[8] - 1
  x[9] - 1
  x[10] - 1
  x[11] - 1
  x[13] - x[14]
```
## Verify claims from the article
To verify Tables 2, 3, and 4 from the article, one can use queries to the database. The example below shows how to generate the line of Table 2 with respect to the matroids of rank 3 on 7 elements. Recall that the uniform matroids are not stored in the database.
The other rows as well as Table 3 can be verified similarly. Note that the database collection for (4,9) is not filled completely yet.
```
julia> t2 = [length([r for r in find(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"3", "data.length_groundset"=>"7","data.dim_s"=>i]))]) for i in ["-1","0","1","2","3","4","5"]]
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
```
julia> t4 = [length([r for r in find(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"4", "data.length_groundset"=>"9","data.dim_r"=>i,"data.dim_s"=>"-1"]))]) for i in ["0","1","2","3","4","5","6"]]
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
```
julia> notequal = find(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"4", "data.length_groundset"=>"9","data.equality_of_realizationspaces"=>"false"]));

julia> length([r for r in notequal])
5399
```
Theorem 4.11 claims that there are at least 5400 matroids of rank 4 on 9 elements with $\mathcal{R}\supsetneq\mathcal{S}$. Since the database does not count the uniform matroid $U_{4,9}$, the claim is verified.

The user can find the selfprojecting matroids for which the computation of the selfprojecting realization space was too costly and did not terminate as follows:
```
julia> notterminated = find(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"3", "data.length_groundset"=>"8","data.dim_s"=>nothing]));

julia> length([r for r in notterminated])
4
```

## Additional Code to verify Tables and Examples

The following explains how to use the code at https://github.com/AlheydisGeiger/selfprojectingGrassmannian to obtain the tables and examples from the paper.

**Tables 2 and 3**

These tables show the distribution of the dimension of the (self-projecting) realization spaces for matroids of rank k on n elements, for which the realization spaces were computed.
The code is included in the file ``generating_tables.jl``. There are two functions, one for the dimensions of the realization spaces, in the article denoted by $\mathcal{R}$, and one for the dimensions of the self-projeting realization spaces $\mathcal{S}$; see Definition 4.8.

The functions need to access the magma output files, stored in the magma directory of this gitHub repository. These files need to be downloaded in order to reproduce the tables.
```
julia> using Oscar
julia> include("your/path/to/generating_tables.jl")
julia> generate_table_content_dimR("your/path/to/rk3on8.out",3,8)
The dimensions of the realization spaces for self-projecting rank 3 matroids on 8 elements are distributed as follows (without the uniform matroid)
[-1   0   1    2    3   4   5   6   7   8]
[ 2   2   5   11   12   9   3   3   1   0]

julia> generate_table_content_dimS("your/path/to/rk3on8.out",3,8)
The dimensions of the self-projecting realization spaces for self-projecting rank 3 matroids on 8 elements are distributed as follows (without the uniform matroid)
[-1   0   1    2    3   4   5   6   7   8]
[ 2   2   5   11   12   9   3   3   1   0]
```

```
julia> include("your/path/to/generating_tables.jl")
julia> generate_table_content_dimR("your/path/to/rk3on8.out",3,8)
The dimensions of the realization spaces for self-projecting rank 3 matroids on 8 elements are distributed as follows (without the uniform matroid)
[-1   0   1    2    3   4   5   6   7   8]
[ 2   2   5   11   12   9   3   3   1   0]

julia> generate_table_content_dimS("your/path/to/rk3on8.out",3,8)
The dimensions of the self-projecting realization spaces for self-projecting rank 3 matroids on 8 elements are distributed as follows (without the uniform matroid)
[-1   0   1    2    3   4   5   6   7   8]
[ 2   2   5   11   12   9   3   3   1   0]
```
When using as inputfiles the magma output files containing only the computations of realization spaces $\mathcal{R}$, the functions need to be endowed with an additional ``_onlyR`` at the end of their names. The inputfiles are similary marked.
```
julia> include("your/path/to/generating_tables.jl")
julia> generate_table_content_dimR_onlyR("your/path/to/opt_3_8_onlyR.out",3,8)
The dimensions of the realization spaces for self-projecting rank 3 matroids on 8 elements are distributed as follows (without the uniform matroid)
[-1   0   1    2    3    4   5   6   7   8]
[ 2   2   5   11   12   11   5   3   1   0]
```


**Table 4**
This table shows the dimensions of realization spaces $\mathcal{R}(M)$ of self-projecting matroids of rank 4 on 9 elements with $\mathcal{S}= \emptyset$. The functions demonstrated below show how to access the identifiers for given $\dim(\mathcal{R})$ and to how to fill the row entries of Table 4.
```
julia> include("your/path/to/generating_tables.jl")
julia> find_all_realizable_not_sp_realizable("realisation_output/rank4/rk4on9_all.out",0)
4-element Vector{Any}:
 5985
 2788
 6265
 7274
julia> length(find_all_realizable_not_sp_realizable("realisation_output/rank4/rk4on9_all.out",5))
124
```


**Example 4.12**

In order to work with the database and/or compute self-projecting realization spaces of matroids in OSCAR, you need to use the developers version of OSCAR on the branch ag/selfprojecting_matroids on <https://github.com/AlheydisGeiger/Oscar.jl/tree/ag/selfprojecting_matroids>.
To reproduce example 4.12 you can access the relevant file from the database.
```
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

Project contributors: Alheydis Geiger, Francesca Zaffalon.

Corresponding author of this page: Alheydis Geiger, 
<a href="mailto:geiger\@mis.mpg.com">geiger\@mis.mpg.de</a>

 
Software used: Magma (V2.27), Julia (Version 1.12.1), OSCAR (version 1.6.0-DEV), 
GNU parallel 20221122


Last updated 10/12/2025.

