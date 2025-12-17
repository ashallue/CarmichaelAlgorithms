# CarmichaelAlgorithms
New project in 2025: algorithms for Carmichael numbers

  1 - Algorithm that answers "Is n a Carmichael number?" in deterministic polynomial time (run-time, but not correctness, depends on ERH).
  
  2 - Algorithm that tabulates all Carmichael numbers less than B asymptotically faster than prior works
  
      2a - The tabulation found here is hard-coded to B = 10^24 and makes several assumptions about the state of inputs
      
      
Problems with prior work (see [the paper](https://link.springer.com/article/10.1007/s40993-024-00598-3) and [the GitHub repo](https://github.com/ashallue/tabulate_car/tree/master)).

  1 - Too many preproducts were created.
  
    1a - The joint work sahred between each process limited the effectiveness of parallelization (see Amdahl's Law). 
    
    1b - What should have been lower-order computations (preproduct generation, cyclic checks) became the bottleneck
    
 2 - The multiple tabulation regimes (one for each count of distinct prime factors) created a complex code-base.
 
    2b - In theory, recursion could have simplified the code base but timing test shows significant slowdown and it was harder to parallelize.
    

There are four main components to be seen here:

  1 - CN_query.cpp is a standalone program that answers "Is n a Carmichael number?"
  
  2 - Preproduct.cpp and Preproduct.h are the main drivers of the tabulation algorithm.
  
    2a - See detailed comments in these that describe their intended usage
    
    2b - It also contains a specialized version of CN_query designed around knowing a significant pre-factored portion of n.
    
  3 - A folder containing tabulation related programs.  
  
    3a - Two major drivers for our tabulation
    
    3b - Miscellany test files and one sample test case that should finish in a reasonable time
    
  4 - A folder containing our modification to Jon Sorenson's bit-packed incremental sieve (see [the paper](https://www.sciencedirect.com/science/article/abs/pii/S002001902400067X) and the [the GitHub repo](https://github.com/sorenson64/soespace))
