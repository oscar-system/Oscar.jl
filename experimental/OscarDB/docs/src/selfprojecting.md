# Realizations of self-projecting matroids
This website accompanies the article "The self-projecting Grassmannian" by Alheydis Geiger and Francesca Zaffalon. It contains the database of (self-projecting) realization spaces of self-projecting matroids of rank k on n elments over characteristic zero for (k,n) in {(2,4),...,(2,12),(3,6),(3,7),(3,8),(4,8),(4,9),(5,10)}.



**Abstract:**  We introduce the self-projecting Grassmannian, an irreducible subvariety of the Grassmannian parametrizing linear subspaces that satisfy a generalized self-duality condition. We study its relation to classical moduli spaces, such as the moduli spaces of pointed curves of genus $g$, as well as to other natural subvarieties of the Grassmannian. We further translate the self-projectivity condition into the combinatorial language of matroids, introducing self-projecting matroids, and we computationally investigate their realization spaces inside the self-projecting Grassmannian.



For the cases {(2,4),(3,6),(4,8),(5,10)} the database stores material from the  article:
Alheydis Geiger, Sachi Hashimoto, Bernd Sturmfels, Raluca Vlad: Self-dual matroids from canonical curves
In: Experimental mathematics, 33 (2024) 4, p. 701-722
DOI: `10.1080/10586458.2023.2239282 <https://dx.doi.org/10.1080/10586458.2023.2239282>`_ ARXIV: https://arxiv.org/abs/2212.05910 CODE: https://github.com/sachihashimoto/self-dual
Project page created: 18/11/2025.

The following explains how to use the code at https://github.com/AlheydisGeiger/selfprojectingGrassmannian to obtain the tables and examples from the paper.


**Tables 2-4**
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

**Example 4.15**



Project contributors: Alheydis Geiger, Francesca Zaffalon.

Corresponding author of this page: Alheydis Geiger, 
<a href="mailto:geiger\@mis.mpg.com">geiger\@mis.mpg.de</a>

 
Software used: Magma (V2.27), Julia (Version 1.12.1), OSCAR (version 0.11.0), GNU Parallel, Macaulay2 (version 1.24.11)

Last updated 18/11/2025.

