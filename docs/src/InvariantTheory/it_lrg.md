```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["it_lrg.md"]
```

# Invariants of Linearly Reductive Groups

Using the notation from the introductory section to this chapter, we now suppose that $\rho: G\to V$ is a *rational* representation
of a *linearly reductive* group $G$ on $V$. 

!!! note
    
    - By the very definition of linear reductivity, there exists a Reynolds operator $\mathcal R: K[V] \to K[V]$. 
    - By Hilbert's celebrated finiteness theorem, $K[V]^G$ is finitely generated as a $K$-algebra.
    - As shown by Hochster and Roberts, $K[V]^G$ is Cohen-Macaulay. 


Omega-process, Derksen's algorithm
