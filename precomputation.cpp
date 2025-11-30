
#include "rollsieve.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <array>

// if the bound or elimination rule is changed drastically
// the use of uint64_t to store preproducts will fail



int main()
{
    uint64_t P = 999'901; // initialize at an odd number
    Rollsieve incremental_sieve( P );
    std::vector< uint64_t > P_fac ;
    
    __int128 B = 1'000'000'000'000;
    B = B * B;
    uint64_t num_jobs =1;
    
    bool is_admissible;
    
    while( P < 1'001'000 )
    {
        is_admissible = false;
        
        incremental_sieve.next();
        incremental_sieve.next();

        P = incremental_sieve.getn();
        incremental_sieve.getlist( P_fac );
        
        if( P_fac.size() == 1 )
        {
            if( P_fac[0] == P && P != incremental_sieve.s )
            {
                is_admissible = true;
            }
        }
        else
        {
            uint64_t sq_free = 1;
            for( auto p : P_fac )
            {
                sq_free *= p;
            }
            // only check on sqaure-free numbers
            if( sq_free == P )
            {
                // set is_admissible to true and see if we should undo that with congruence checks
                is_admissible = true;
                std::sort( P_fac.begin(), P_fac.end() );
                uint16_t p_count = P_fac.size();
                uint16_t i = 0;
                uint16_t j = 1;
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
        
        if( is_admissible )
        {
            // generally admissible but now we need to make sure it is bounds admissible
            // that is, Pp^3 < B
            
            bool bounds_admiss = ( (__int128) P  <= ( B / ( (__int128) P_fac.back() * (__int128) P_fac.back()* (__int128) P_fac.back() ) ) );
             
            if( bounds_admiss )
            {
                num_jobs++;
                
                std::cout << P << " =  " ;
                uint64_t L = 1;
                for( auto p : P_fac )
                {
                    std::cout << p << " ";
                    L = std::lcm( L , p - 1 );
                }
                std::cout << ", " << L << ", " << std::max( 125'000'000/P, P_fac.back()) <<  std::endl;
                 
            }
        }
    }
    std::cout << " You have " << num_jobs << " admissible preproducts to consider." << std::endl;
}
