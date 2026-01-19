#include "Preproduct.h"
#include "IncrementalSieve/rollsieve.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdio.h>
#include <array>
#include <vector>
#include <chrono>

using namespace std::chrono;

/*
This file is an independent program designed to do smaller serial tabulations for testing 
and timing purposes.
*/

// compiled and run with:
// mpic++ -o parallel parallel_large_P.cpp Preproduct.o rollsieve.o -lgmp -O3
// mpirun --hostfile hostfile.txt -n 156 ./parallel &

int main(int argc, char * argv[])
{
    int my_rank = 0;            // my CPU number for this process
    int proc = 1;               // number of CPUs that we have

    int B_pow = 17;             // B will be 10 to this power
    
    //std::string file_name = "nbd4_large_" + std::to_string(my_rank) + ".txt";
    std::string file_name = "cars_large10to" + std::to_string(B_pow) + ".txt";
    //const int node_num = 4;
    
    uint64_t P = 3; // initialize at an odd number
    Rollsieve incremental_sieve( P );
    std::vector< uint64_t > P_fac ;
    
    __int128 B = pow(10, B_pow);
    uint64_t X;
    double one_third = 1.0 / 3;
    X = ceil(pow(B, one_third));
    
    uint64_t num_jobs =1;
    
    bool is_admissible;
    uint64_t count_admissible = 0;
    
    Preproduct large_P = Preproduct();

    auto start_large = high_resolution_clock::now();
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
                // check admissibility
                // we could check gcd( P, phi(P) ) == 1
                // but that requires computing phi(P) and then the Euclidean algorithm
                // instead, we check admissibility by congruences with the primes
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
                //num_jobs++;
                //if( ( my_rank + ( node_num * 112 ) ) == ( num_jobs % 560 ) )
                //{
                 //   std::cout << P ;
                    count_admissible++;
                
                    uint64_t L = 1;
                    for( auto p : P_fac )
                    {
                        // std::cout << p << " ";
                        L = std::lcm( L , p - 1 );
                    }
                    large_P.initializing( P, L, std::max( X/P, P_fac.back() ) );
                    large_P.CN_multiples_of_P( file_name );
                //}
            }
        }
    }

    auto end_large = high_resolution_clock::now();

    auto duration_large = duration_cast<seconds>(end_large - start_large);
    std::cout << "Timing for large preproduct case with new code base, is: " << duration_large.count() << "\n\n";
    std::cout << "Count of bounds admissible preproducts: " << count_admissible << "\n";
    return 0;
}
