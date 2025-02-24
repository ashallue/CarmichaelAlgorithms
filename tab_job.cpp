#include "Preproduct.h"

#include <iostream>
#include <algorithm>
#include <gmp.h>
#include <numeric>
#include <chrono>

int main(void) {
    
    Preproduct P0;
    // P0.initializing( 599266767, 890750, 991 );
    P0.initializing( 6682828353, 2289560, 13 );
    
    std::cout << "LP for this is " << sizeof( unsigned long int ) << std::endl;
    
    std::cout << "Initializing P : " ;
    gmp_printf ("%Zd = ", P0.P );
    
    for( int i = 0; i < P0.P_primes.size(); i++)
    {
        std::cout << P0.P_primes[i] << " * "  ;       
    }
    std::cout << P0.P_primes[P0.P_primes.size() - 1 ] << std::endl ;
    
    std::cout << "Initializing Lambda : " ;
    gmp_printf ("%Zd = ", P0.L );
    std::cout << std::endl;
   
    /*  A good test case
    Preproduct P1;
    P1.initializing( 515410417841, 115920, 50 );
    P1.CN_search();
    */
     
    // the triple 132582235, 91872, 941 is in the working jobs with the new rule
    // this means we have between 2 and 6 primes to append to find a CN
    // for all primes in ( 941, (B/P)^(1/3) ) = (941, 196112) we append a prime and do CN_search
    // for all primes in ( 196112, (B/P)^(1/2) ) = ( 196112, 86847502) we append a prime and look for exactly one prime
    // this represents the first two CN_xearches invoked by the above:
    

    auto start1 = std::chrono::high_resolution_clock::now();

   
    uint64_t P = 132582235;
    uint64_t L = 91872;
    
    uint64_t P1;
    uint64_t L1;
    
    uint64_t p64 = 947;
    mpz_t p ;
    mpz_init_set_ui( p, p64 );
    
    Preproduct PPP;
    PPP.initializing( 132582235, 91872, 900  );
   
    // this example is *not* ideal, we refactor P by trial division every time
    // we need to have p2-1 in factored form in order to use appending call
    // as you can see below, for ease, I'm just using gmp's next prime
    // and am just re-initializing each time.
    Preproduct P_testing;

    while( p64 < 196112 )
    {
        mpz_nextprime( p, p);
        p64 = mpz_get_ui( p );
        if( PPP.is_admissible_modchecks( p64 ) )
        {
            //std::cout << "I found an admissible prime " << p64 << std::endl;
            P1 = P*p64;
            L1 = L*((p64-1)/std::gcd( L , p64-1 ) );
            P_testing.initializing ( P1, L1, p64 );
            P_testing.CN_search( "cars.txt" );
            
        }
    }
    mpz_clear( p );

    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

    std::cout << "Time taken: " << duration1.count() << " microseconds" << std::endl;
    
    /*
    auto start2 = std::chrono::high_resolution_clock::now();
    
    
    Rollsieve r(941);
    p64 = r.getn();
    std::vector< uint64_t > factors;
    
    Preproduct PP;
    PP.initializing( 132582235, 91872, 900  );
       
    while( p64 < 196112 )
    {
        if( r.isnextprime() && PP.is_admissible_modchecks( p64+1 ) )
        {
            r.getlist( factors );
            std::sort( factors.begin(), factors.end() );

            P_testing.appending( &PP, p64+1, factors );
            P_testing.CN_search( "cars.txt" );
        }
        r.next();
        p64 = r.getn();
    }
    
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);

    std::cout << "Time taken: " << duration2.count() << " microseconds" << std::endl;

    /*
     This run completes in 4 minutes and 30 seconds on thomas.butler.edu and finds these:
     433493717815335774335905 =  5 17 23 73 929 2377 7129 10267 18793
     72425097332690148535105 =  5 17 23 73 929 5347 577 673 263089
     433493717815335774335905 =  5 17 23 73 929 7129 2377 10267 18793
     433493717815335774335905 =  5 17 23 73 929 10267 2377 7129 18793
     433493717815335774335905 =  5 17 23 73 929 18793 2377 7129 10267
     725906640592907462305 =  5 17 23 73 929 21577 643 394633
     
     The fact that these are found duplicated only indicates our failure to implement the append_bound-related checks in CN_factorization
     */
   
    
    
  
    
    // P0.CN_search(1873371784);
    //P0.CN_search(149637241475922);

    
    return 0;
}

