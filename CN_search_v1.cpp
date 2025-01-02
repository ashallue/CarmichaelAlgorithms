// compiled with  g++ CN_search.cpp -lgmp -O3

// no wheel
// no sieving


#include <gmp.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>
#include <boost/dynamic_bitset.hpp>

int main()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();
      
    // set search bound, set as 10^24
    mpz_t bound;
    mpz_init_set_ui( bound, 10 );
    mpz_pow_ui( bound, bound, 24 );

    // set preproduct
    mpz_t P;
    mpz_init_set_ui( P, 515410417841 );

    // set L = lambda(P)
    mpz_t L;
    mpz_init_set_ui( L, 115920 );
    
    // We do not need other parts of the preproduct class for this approach
    
    // compute r^* = P^{-1} mod L
    // this is the start of  R = (r^* + kL) w/ k = 0
    mpz_t r_star;
    mpz_init(r_star);
    mpz_invert(r_star, P, L);

    mpz_t base2;
    mpz_t base3;
    mpz_init_set_ui( base2, 2 );
    mpz_init_set_ui( base3, 3 );
    
    mpz_t fermat_result;
    mpz_init( fermat_result );
  
    // This is the start of n = Pr^* + kPL w/ k = 0
    mpz_t n;
    mpz_init(n);
    mpz_mul( n, P, r_star);

    // n is an arithemtic progression with common different PL
    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );
    
    while( mpz_cmp( n , bound ) < 0 )
    {
        mpz_powm( fermat_result,  base2,  n, n); // 2^n mod n
            if( mpz_cmp( fermat_result, base2 ) == 0 )  // check if 2 = 2^n mod n
            {
                mpz_powm( fermat_result,  base3,  n, n); // 3^n mod n
                if( mpz_cmp( fermat_result, base3 ) == 0 )  // check if 3 = 3^n mod n
                {
                    // check CN here
                    std::cout << "n = " ;
                    gmp_printf( "%Zd", n);
                    std::cout << " is a base-2 and base-3 Fermat psp." << std::endl;
                }
            }
        // next element in the arithmetic progression
        mpz_add( n, n, PL);
    }

    mpz_clear( P );
    mpz_clear( L );
    mpz_clear( r_star );
    mpz_clear( n );
    mpz_clear( PL );
    mpz_clear( base2 );
    mpz_clear( base3 );
    mpz_clear( fermat_result );
    mpz_clear( bound );

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count() << "ms\n";
}
