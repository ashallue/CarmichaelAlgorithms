# Andrew's Notes

This file records testing, timings, and other things I want to remember like scripting commands.

# Small Tabulations and Timings

While the project is designed for B = 10^24 only, there is of course a need to do smaller tabulations for timings and other checks.  The file for this is small_tabulation.cpp, which compiles to the executable small_tab.  This file combines work found in precomputation.cpp and Preproduct.cpp into a single serial computation.  These runs were performed with the following steps.  

* Set B_pow.  The upper bound B will be 10 to this power.  Note X is automatically set to B^{1/3}.
* In Preproduct.cpp, make sure the #define TEST preprocessing flag is set
* In Preproduct.cpp constructor, set the testing bound to the same value as B_pow
* In Preproduct::CN_multiples_of_P, set the testing value for X to B^{1/3}
* In that same function, just above definition of X, you can choose a flag for early abort or not.  [More details needed]

Now compile and run small_tab.  To do comparisons of files, it is a good idea to sort and remove duplicates.  For example, 

sort -u -n -k1 cars_large10to15.txt > cars15_sorted.txt