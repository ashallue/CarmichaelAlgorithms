//  This file was responsible for tabulating n = PR with omega(P) > 2 and P > 91548
//  It uses the partitioning strategy of Subsection 5.2 of https://arxiv.org/html/2506.09903v3
//      Comments interspersed below show how the triples are computed
//  It was ran on 5 112-core processors and had a separate compiler per machine
//  lines 33-34 changed depended on depending on the node: {0, 1, 2, 3, 4}.  
//  Hard-coded constants for 112 and 5*112 = 560 are found below.

#include "Preproduct.h"
#include "rollsieve.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdio.h>
#include <array>
#include <vector>
#include "mpi.h"

// compiled and run with:
// mpic++ -o parallel parallel_large_P.cpp Preproduct.o rollsieve.o -lgmp -O3
// mpirun -n 112 ./parallel &

int main(int argc, char * argv[])
{
    int my_rank;            // my CPU number for this process
    int proc;               // number of CPUs that we have

    MPI_Init(&argc, &argv);                     // Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    // Find my rank
    MPI_Comm_size(MPI_COMM_WORLD, &proc);       // Find out the number of processes!

    std::string file_name = "nbd4_large_" + std::to_string(my_rank) + ".txt";
    const int node_num = 4;

    // initialize the starting preproduct to an odd number:
    uint64_t P = 91549; 
    Rollsieve incremental_sieve( P );
    std::vector< uint64_t > P_fac ;

    // the constants associted to our two bounds
    __int128 B = 1'000'000'000'000;
    B = B * B;
    const uint64_t X = 125'000'000;

    // counter for parallel distribution
    uint64_t num_jobs = 1;
    
    bool is_admissible;
    
    Preproduct large_P = Preproduct();
    
    while( P < X )
    {
        is_admissible = false;
        
        // increment to next odd number:
        incremental_sieve.next();
        incremental_sieve.next();

        P = incremental_sieve.getn();
        incremental_sieve.getlist( P_fac );
        
        // checks that P is a prime and not a prime power
        if( P_fac.size() == 1 )
        {
            if( P_fac[0] == P && P != incremental_sieve.s )
            {
                is_admissible = true;
            }
        }
        // so the number is the product of at least two distinct primes
        else
        {
            uint64_t sq_free = 1;
            for( auto p : P_fac )
            {
                sq_free *= p;
            }
            // check that the number is square free
            if( sq_free == P )
            {
                // set is_admissible to true and see if we should undo that with congruence checks
                is_admissible = true;
                std::sort( P_fac.begin(), P_fac.end() );
                uint16_t p_count = P_fac.size();
                uint16_t i = 0;
                uint16_t j = 1;
                // check admissibility with congruences
                while( is_admissible && i < (p_count - 1) )
                {
                    while( is_admissible && j < p_count )
                    {
                        is_admissible = ( is_admissible && ( 1 != P_fac[j] % P_fac[i] ) );
                        j++;
                    }
                    i++;
                    j = i+1;
                }
            }
        }
        
        // entering this branch means that P is cyclic
        if( is_admissible )
        {
            // check bounds admissible:  Pp^3 < B
            bool bounds_admiss = ( (__int128) P  <= ( B / ( (__int128) P_fac.back() * (__int128) P_fac.back()* (__int128) P_fac.back() ) ) );
             
            if( bounds_admiss )
            {
                // we parallelize by count of bounds admissible preproducts
                num_jobs++;
                if( ( my_rank + ( node_num * 112 ) ) == ( num_jobs % 560 ) )
                {
                    // compute lambda(P)
                    uint64_t L = 1;
                    for( auto p : P_fac )
                    {
                        L = std::lcm( L , p - 1 );
                    }
                    // We now have the tuple ( P, lambda(P), max{ X/P, p } )
                    large_P.initializing( P, L, std::max( X/P, P_fac.back() ) );
                    large_P.CN_multiples_of_P( file_name );
                }
            }
        }
    }

    MPI_Finalize();

    return 0;
}

/*
 
May 29 10:00am:   ~93k starting point
May 30 10:00am:  ~106k (avg of recent 10)
May 31 10:00am:  ~119k
June 2 10:00am:  ~134k
June 3 10:00am:  ~154k
June 4 10:00am:  ~182k
June 5 10:00am:  ~193k (about 2.07x from start in 1 week)
June 6 10:00am:  ~224k
June 9 10:00am:  ~292k
June 10 10:00am: ~322k
June 11 10:00am: ~345k
June 12 10:00am: ~388k (about 2.01x from last week)
June 13 10:00am: ~438k
June 16 10:00am: ~555k
June 17 10:00am: ~688k
June 18 10:00am: ~740k
June 19 10:00am: ~814k (about 2.10x from last week)
June 20 10:00am: ~866k
June 23 10:00am: ~1.25m
June 24 10:00am: ~1.36m
June 25 10:00am: ~1.56m
June 26 10:00am: ~1.78m (about 2.2x from last week)
...
July 7  10:00am: ~7.49m (averaging 20 now)
July 8  10:00am: ~7.84m
July 9  10:00am: ~11.5m
July 10 10:00am: ~12.4m (about 7x from two weeks ago ~ 2.64^2)
July 11 10:00am: ~15.9m
July 12 10:00am: ~21.1m
July 14 10:00am: ~30.3m (1 of 560 finished at 4:52 am)
July 15 10:00am: ~39.8m (4 of 560 finished)
July 16 10:00am: ~49.9m (20 of 560 finished)
July 17 10:00am: ~71.5m (69 of 560 finished)
July 18 10:00am: 142 of 560 finished
July 19 10:00am: 239 of 560 finished
July 20 10:00am: 366 of 560 finished
July 21 10:00am: 451 of 560 finished
July 22 10:00am: 511 of 560 finished
July 23 10:00am: 542 of 560 finished
July 24 10:00am: 553 of 560 finished
July 25 10:00am: 555 of 560 finished
 
 */
