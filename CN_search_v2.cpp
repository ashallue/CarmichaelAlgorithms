// compiled with  g++ CN_search.cpp -lgmp -O3

// yes wheel
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
    const uint64_t L_len = 5;
    uint64_t L_distinct_primes[L_len] = { 2, 3, 5, 7, 23 };
    
    // controls the size of the wheel
    const uint32_t cache_bound = 150'000;

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

    
    // This is the start of n = Pr^* + kPL w/ k = 0
    mpz_t n;
    mpz_init(n);

    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );

    // position i has the truth value of the statement "(2i + 3) is prime"
    std::bitset<256> small_primes{"0010010100010100010000010110100010010100110000110000100010100010010100100100010010110000100100000010110100000010000110100110010010010000110010110100000100000110110000110010010100100110000110010100000010110110100010010100110100110010010110100110010110110111"};
    
    // turn off the bits corresponding to primes dividing L
    // this will need to be modified so that we don't try to turn off bits out of range
    // that is, we need to be carefuly if L has a prime dividing it that exceeds 513
    for( uint16_t i = 1; i < L_len; i++ )
    {
        small_primes [ ( L_distinct_primes[i] - 3) /2 ] = 0;
    }

    // we need to append primes to L so that B/(PL*appended_primes) < cache_bound
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    mpz_cdiv_q( cmp_bound, bound, PL);
    
    // we need to store the primes that we append
    std::vector< uint16_t > primes_lifting_L;
    
    // L_lift will store the product of primes
    // and also tells the mximum number of spokes on the wheel
    uint64_t L_lift= 1;
    uint16_t prime_index = 0;
    
    while( mpz_cmp_ui( cmp_bound, cache_bound ) > 0  )
    {
        if( small_primes[ prime_index ] )
        {
            uint16_t p = 2*prime_index + 3;
            L_lift *= p;
            primes_lifting_L.push_back( p );
            mpz_mul_ui( PL, PL, p );
            mpz_cdiv_q( cmp_bound, bound, PL);
        }
        prime_index++;
    }
  
    std::cout << "The number of residues per spoke is " <<
    gmp_printf( "%Zd", cmp_bound );
    std::cout << std::endl;
    std::cout << "We are using the primes " ;
    for( auto p : primes_lifting_L )
    {
        std::cout << p << " " ;
    }
    std::cout << "in our wheel" << std::endl;
    
    // This loop is the outer wheel
    for( uint64_t m = 0; m < L_lift; m ++)
    {
        bool enter_loop = true;
        // checks to make sure that r_star + m*L is not divisible by lifted_primes
        for( auto p : primes_lifting_L )
        {
            enter_loop = ( enter_loop && ( mpz_divisible_ui_p( r_star, p ) == 0 ) );
        }
        
        // this represent a spoke not divisibly by a wheel prime
        if( enter_loop )
        {
            // this is where n is initialized to check candidates on a spoke of a wheel
            mpz_mul( n, P, r_star);
            while( mpz_cmp( n , bound ) < 0 )
            {
                mpz_powm( fermat_result, base2, n, n); // 2^n mod n
                if( mpz_cmp( fermat_result, base2 ) == 0 )  // check if 2 = 2^n mod n
                {
                    mpz_powm( fermat_result, base3, n, n); // 3^n mod n
                    if( mpz_cmp( fermat_result, base3 ) == 0 )  // check if 3 = 3^n mod n
                    {
                        // CN check here.
                        std::cout << "n = " ;
                        gmp_printf( "%Zd", n);
                        std::cout << " is a base-2 and base-3 Fermat psp." << std::endl;
                    }
                }
                // next element in the arithmetic progression
                mpz_add( n, n, PL);
            }
        }
        
        // we do L_lift steps in arithmetic progression
        // each step represents one spoke" on a wheel
        mpz_add( r_star, r_star, L);
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
