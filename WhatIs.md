# Oscar - what is the idea

The ultimate goal for Oscar should be to compete with (and ideally beat) Magma and Sage in our areas of
expertise.

Oscar is the innovative, next generation Computer Algebra System

Oscar is written in Julia, but is not Julia, nor can it be. 
 Some examples to illustrate what that means:
  Julias matrices are (C) arrays (of arbitrary dimension), parametrized by
  the type of the entries (apart from banded, sparse, ... special matrices).
  In the numerical world, the type mostly defines the representation of an object
  - double and variations
  - complex
  - BigFloat
  - Int
  - BigInt

  In algebra, this is either not true or terribly inefficient (or impossible)
  Take
    Z/nZ integers modulo n, and matrices over it
   
   Then either
    n is part of the type -> every function is recompiled for every n - which kills
    all modular (Chinese remainder theorem  (CRT) based) algorithms

    n is not part of the type, then it need to be elsewhere, ie. in the parent,
    every element, additional arguments, ...

    Furthermore, if n is BigInt (fmpz), so no bittype, then it cannot be part of the type.

  for non-empty matrices, this does not matter as the entries have enough information. For
  empty matrices....

  So: offshot: we simply cannot use normal Julia infrastructure in many places, so please
  use ours. If functions are missing then
  - add them
  - or tell us

The main point about Oscar is for non-experts to allow access, thus it needs to
work as you would expect as the youngest student knowing the objects. Generically, 
Oscar is supposed to follow mathematical conventions. In case of several applying,
the more general one. Oscar is not just supposed to be for experts in a small area
but should support a wide range, hence expect that typical mathematical shortcuts
do not work.
As an example: in classical number theory there is the convention that many definitions
that are trivial for fields are silently applied to the ring of integers. One
speaks of the unit group of the number field, meaning the unit group of the ring
of integers. On paper, the meaning of a term or object is defined by the context. In 
Oscar, this does not work: the meaning has to be part of either
 - the object
 - or the question

  Let alpha be an element explicitly constructed as an element of the number field,
  not of the ring of integers, then

   - is_unit will just test if it is non-zero (unit in a field as a special type of ring)
   - is_unit_in_ring_of_integers would supply the context for the other interpretation.

  In Oscar, this context is mostly supplied by the type of the object and possibly the
  parent, e.g. the ring/ field/ group we're in.

