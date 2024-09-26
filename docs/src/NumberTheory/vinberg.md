```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Vinberg's algorithm

 A $\textit{Lorentzian lattice} L$ is an integral $\Z$-lattice of signature $(s_+, s_-)$ with $s_+=1$ and $s_->0$. 
 A $\textit{root}$ of $r \in L$ is a primitive vector s.t. reflection in the hyperplane $r^\perp$ maps $L$ to itself.
 $L$ is called $\textit{reflective} if the set of fundamental roots $\~{R}(L)$ is finite.

 See for example [Tur18](@cite) for the theory of Arithmetic Reflection Groups and Reflective Lorentzian Lattices.

## Description 
 The algorithm constructs a fundamental polyhedron $P$ for a Lorentzian lattice $L$ by computing its $\textit{fundamental roots} $r$, i.e. the roots $r$ which are perpendicular to the faces of $P$ and which have inner product at least 0 with the elements of $P$.
 Choose $v_0$ in $L$ primitive with $v_0^2 > 0$ as a point that $P$ should contain.

 Let $Q$ be the corresponding Gram matrix of $L$ with standard basis. A vector $r$ is a fundamental root, if
 - the vector $r$ is primitive,
 - reflection by $r$ preserves the lattice, i.e. $\frac{2}{r^2}*r*Q$ is an integer matrix,
 - the pair $(r, v_0)$ is positive oriented, i.e. $(r, v_0) > 0$,
 - the product $(r, \~{r}) \geq \ 0$ for all roots $\~{r}$ already found.
 This implies that $r^2$ divides $2*i$ for $i$ being the level of $Q$, i.e. the last invariant of the Smith normal form of $Q$. 

 $P$ can be constructed by solving $(r, v_0) = n$ and $r^2 = k$ by increasing order of the value $\frac{n^2}{k}$ and $r$ satisfying the above conditions.

 If $v_0$ lies on a root hyperplane, then $P$ is not uniquely determined.
 In that case we need a direction vector $v_1$ which satisfies $(\~{v}, v_1) \neq 0$ 
 for all possible roots $\~{v}$ with $(v_0, \~{v}) = 0$  

 With $v_0$ and $v_1$ fixed $P$ is uniquely determined for any choice of root lengths and maximal distance $(v_0, r)$.
 We choose the first roots $r$ by increasing order of the value $\frac{(\~{r}, v_1)}{r^2}$ for all possible roots $\~{v}$ with $(v_0, \~{v}) = 0$.
 For any other root length we continue as stated above.
 
 For proofs of the statements above and further explanations see [Vin75](@cite).

 ## Function
 
 ```@docs
 vinberg_algorithm(Q::ZZMatrix, upper_bound::ZZRingElem)
 ```


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en).
* [Stevell Muller](https://www.math.uni-sb.de/ag/brandhorst/index.php?option=com_content&view=article&id=30:muller&catid=10&lang=de&Itemid=104),

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).