# CarmichaelAlgorithms
New project in 2024: algorithms for tabulating Carmichael numbers, among others.
Joint with Jonathan Webster.

One issue with our 2024 ANTS paper (and the associated codebase github.com/ashallue/tabulate_car) is that for the large preproduct case, 2-3 gcd computations were required at every entry into a nested loop.  We now believe that the generation of preproducts was the bottleneck in that code, and have new ideas to remove those gcd computations, hopefully resulting in significantly faster tabulation of this case.  Having received critique of the previous codebase, we decided to start fresh for this project.
