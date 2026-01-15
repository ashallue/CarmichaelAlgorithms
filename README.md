# CarmichaelAlgorithms
New project in 2024: algorithms for tabulating Carmichael numbers, among others.
Joint with Jonathan Webster!

One issue with our 2024 ANTS paper (and the associated codebase github.com/ashallue/tabulate_car) is that for the large preproduct case, 2-3 gcd computations were required at every entry into a nested loop.  We now believe that the generation of preproducts was the bottleneck in that code, and have new ideas to remove those gcd computations, hopefully resulting in significantly faster tabulation of this case.  Having received critique of the previous codebase, we decided to start fresh for this project.

Another issue with the ANTS 2024 tabulation was the amount of duplicated work.  The way we parallelized had each processer duplicate work.  Tabulating by the count of distinct prime divisors duplicates work, too.  E.g., many of the admissibility checks done for the d=6 case are just repeated again for the d=7 case. 

I'm going to write some things that really shouldn't be in a readme.  We will need a makefile.  We will need classes.  What should we call those classes?

On the topic of ignoring files, see https://www.atlassian.com/git/tutorials/saving-changes/gitignore
A potential problem: this is in the repository, so I think it would be ignored for all collaborators.  I don't yet know how to get a single collaborator to ignore a file that is tracked by a different collaborator.

Note that this project is designed only for a 10^24 tabulation.  This includes hard-coded bounds (upper and cross-over), and integer size choices that only make sense if the upper bound is 10^24.

# Small Tabulations and Timings

While the project is designed for B = 10^24 only, there is of course a need to do smaller tabulations for timings and other checks.  The file for this is small_tabulation.cpp, which compiles to the executable small_tab.  This file combines work found in precomputation.cpp and Preproduct.cpp into a single serial computation.  To reproduce these check runs, here's what needs to be done:

* Set B_pow.  The upper bound B will be 10 to this power.  Note X is automatically set to B^{1/3}.
* In Preproduct.cpp, make sure the #define TEST preprocessing flag is set
* In Preproduct.cpp constructor, set the testing bound to the same value as B_pow
* In Preproduct::CN_multiples_of_P, set the testing value for X to B^{1/3}
* In that same function, just above definition of X, you can choose a flag for early abort or not.  [More details needed]

Now compile and run small_tab.  To do comparisons of files, it is a good idea to sort and remove duplicates.  For example, 

sort -u -n -k1 cars_large10to15.txt > cars15_sorted.txt