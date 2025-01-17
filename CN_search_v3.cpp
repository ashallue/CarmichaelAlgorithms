// compiled with  g++ CN_search.cpp -lgmp -O3

// no wheel
// yes sieving

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
    const uint64_t L_len = 5;
    uint64_t L_distinct_primes[L_len] = { 2, 3, 5, 7, 23 };

    // all primes less than append_bound are used in sieving
    uint32_t append_bound = 100;
    
    // compute r^* = p^{-1} mod L
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

    // n will be used a temp variable for a bit
    mpz_t n;
    mpz_init(n);

    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );
    
    // position i has the truth value of the statement "(2i + 3) is prime"
    std::bitset<256> small_primes{"0010010100010100010000010110100010010100110000110000100010100010010100100100010010110000100100000010110100000010000110100110010010010000110010110100000100000110110000110010010100100110000110010100000010110110100010010100110100110010010110100110010110110111"};
    // fix append_bound to be consistent with the bitset
    append_bound = (append_bound - 3)/2;
        
    // turn off the bits corresponding to primes dividing L
    // this will need to be modified so that we don't try to turn off bits out of range
    // that is, we need to be carefuly if L has a prime dividing it that exceeds 513
    for( uint16_t i = 1; i < L_len; i++ )
    {
        small_primes [ ( L_distinct_primes[i] - 3) /2 ] = 0;
    }

    // create bit-vector for sieving here:
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    mpz_cdiv_q( cmp_bound, bound, PL);
    uint64_t cmp_bound64 = mpz_get_ui( cmp_bound );
    boost::dynamic_bitset<> spoke_sieve( cmp_bound64 );
    spoke_sieve.reset();
    
    // now sieve
    // the "10" needs to be related to append_bound found in preproduct.h
    mpz_t small_prime;
    mpz_init( small_prime );
    uint16_t prime_index = 0;
    std::cout << "we sieved by " ;
    while( prime_index < append_bound )
    {
        // r_star + k*L = 0 mod p
        // implies k = -r*L^{-1} mod p
       if( small_primes[ prime_index ] )
        {
            int p = 2*prime_index + 3;
            std::cout << p << " " ;
            mpz_set_ui( small_prime, p );
            mpz_invert( n, L, small_prime );  // n has L^{-1} mod p
            mpz_neg( n, n);  // n has -(L)^{-1} mod p
            mpz_mul( n, n, r_star ); // muliply by r_star
            mpz_mod( n, n, small_prime );  // reduce modulo p
            uint64_t k = mpz_get_ui( n );
            
            // we just found the starting point for k, now sieve:
            while( k < cmp_bound64 )
            {
                spoke_sieve[k] = 1;
                k += p;
            }
        }
        prime_index++;
    }
    std::cout << std::endl;
    mpz_clear( small_prime );
          
    // sieve is set up
    // do normal loop but only enter on unmarked items
    // compare to v1
    mpz_mul( n, P, r_star);
    uint32_t k = 0;
    while( mpz_cmp( n , bound ) < 0 )
    {
        if( spoke_sieve[k] == 0 )
        {
            mpz_powm( fermat_result,  base2,  n, n); // 2^n mod n
            if( mpz_cmp( fermat_result, base2 ) == 0 )  // check if 2 = 2^n mod n
            {
                mpz_powm( fermat_result,  base3,  n, n); // 3^n mod n
                if( mpz_cmp( fermat_result, base3 ) == 0 )  // check if 3 = 3^n mod n
                {
                    std::cout << "n = " ;
                    gmp_printf( "%Zd", n);
                    std::cout << " is a base-2 and base-3 Fermat psp." << std::endl;
                }
            }
        }
        k++;
        mpz_add( n, n, PL);
    }
    
          
    mpz_clear( cmp_bound );
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
