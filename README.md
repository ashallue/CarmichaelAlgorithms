# CarmichaelAlgorithms
New project in 2025: algorithms for Carmichael numbers:
  - An algorithm that answers "Is n a Carmichael number?" in deterministic polynomial time (run-time, but not correctness, depends on ERH).
  - An algorithm that tabulates all Carmichael numbers less than B asymptotically faster than prior works.
    - The tabulation found here is hard-coded to B = 10^24 and makes several assumptions about the state of inputs.
      
      
The prior work had some practical problems (see [the paper](https://link.springer.com/article/10.1007/s40993-024-00598-3) and [the GitHub repo](https://github.com/ashallue/tabulate_car/tree/master)):
  - Too many preproducts were created.
    - The joint work shared between each process limited the effectiveness of parallelization (see Amdahl's Law). 
    - What should have been lower-order computations (preproduct generation, cyclic checks, lambda(P) computations, etc.) became the bottleneck.
  - The multiple tabulation regimes (one for each count of distinct prime factors) created a complex code-base.
    - In theory, recursion simplifies the code base but timing tests showed significant slowdown compared to iterative approach.

There are four main components in this repository:
  - CN_query.cpp is the standalone program that answers "Is n a Carmichael number?"
  - Preproduct.cpp and Preproduct.h are the main drivers of the tabulation algorithm.
    - See detailed comments in these that describe their intended usage.
    - It also contains a specialized version of CN_query designed around knowing a significant pre-factored portion of n.
  - A folder containing tabulation related programs.  
    - Two major drivers for our tabulation.  These created the necessary inputs for the tabulation.  
    - Miscellany test files and one sample test case that should finish in a reasonable time.
  - A folder containing our modification to Jon Sorenson's bit-packed incremental sieve of Eratosthenes (see [his paper](https://www.sciencedirect.com/science/article/abs/pii/S002001902400067X) and his [GitHub repo](https://github.com/sorenson64/soespace)).
