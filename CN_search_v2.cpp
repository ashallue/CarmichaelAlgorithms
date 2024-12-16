// compiled with  g++ CN_search.cpp -lgmp -O3

#include <gmp.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>

int main()
{
    // the below has the following quantities hard-coded:
    // The upper bound of the computation:      lines 33,34
    // The quantity P:                          line 38
    // The quantity L = lambda(P)               lines 44,45
    // The array holding the primes dividing L  lines 46
    // e = v_2( LCM( P-1, L0 ) )                line 54
    
    // version 1: 401 s, mpz_class version
    // version 2: 376 s, machine words where possible, import to mpz_t for exponentiation
    // version 3: 338 s, all arithmetic in mpz_t
    // so we go with all arithmetic in mpz_t
    
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
    uint64_t L_distinct_primes[5] = { 2, 3, 5, 7, 23 };
    
    // compute r^* = p^{-1} mod L
    // this is the start of  R = (r^* + kL) w/ k = 0
    mpz_t r_star;
    mpz_init(r_star);
    mpz_invert(r_star, P, L);

    // will need more bases later
    // use bases from the prime divisors of L
    mpz_t base2;
    mpz_init_set_ui( base2, 2 );
    mpz_t base3;
    mpz_init_set_ui( base3, 3 );
    
    mpz_t base_2_fermat;
    mpz_init( base_2_fermat );
    mpz_t base_3_fermat;
    mpz_init( base_3_fermat );

    
    // This is the start of n = Pr^* + kPL w/ k = 0
    mpz_t n;
    mpz_init(n);

    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );

    // Here:  construct a new PL so that B/PL < 10^5
    // We hard-code it for this example:

    // { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    // pull from the small primes in order:
    uint32_t small_sieve_primes[4] = { 11, 13, 17 };
    
    // lifting by 11, 13, 17
    mpz_mul_ui( PL, PL, small_sieve_primes[0] );
    mpz_mul_ui( PL, PL, small_sieve_primes[1] );
    mpz_mul_ui( PL, PL, small_sieve_primes[2] );
    
    // Now B/PL < 10^5
    // Now we can sieve bit arrays of maximal length 10^5 bits for primes not dividing PL/P
    
    for( uint64_t m = 0; m < 11*13*17; m ++)
    {
        if( mpz_divisible_ui_p( r_star, small_sieve_primes[0] ) == 0
            && mpz_divisible_ui_p( r_star, small_sieve_primes[1] ) == 0
            && mpz_divisible_ui_p( r_star, small_sieve_primes[2] ) == 0
           )
        {
            mpz_mul( n, P, r_star);
            
            while( mpz_cmp( n , bound ) < 0 )
            {
                mpz_powm( base_2_fermat,  base2,  n, n); // 2^n mod n
                if( mpz_cmp( base_2_fermat, base2 ) == 0 )
                {
                    mpz_powm( base_3_fermat,  base3,  n, n); // 3^n mod n
                    if( mpz_cmp( base_3_fermat, base3 ) == 0 )
                    {
                        // n is now base 2 and a base 3 Fermat psp
                        // invoke CN factorization algorithm here
                        std::cout << "n = " ;
                        gmp_printf( "%Zd", n);
                        std::cout << " is a base-2 and base-3 Fermat psp." << std::endl;
                    }
                }
                mpz_add( n, n, PL);
            }
        }
        
        mpz_add( r_star, r_star, L);  //next element in arithmetic progression
    }
        
    mpz_clear( P );
    mpz_clear( L );
    mpz_clear( r_star );
    mpz_clear( n );
    mpz_clear( PL );
    mpz_clear( base2 );
    mpz_clear( base3 );
    mpz_clear( base_2_fermat );
    mpz_clear( base_3_fermat );
    mpz_clear( bound );

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count() << "ms\n";



}
