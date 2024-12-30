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
    const uint64_t L_len = 5;
    uint64_t L_distinct_primes[L_len] = { 2, 3, 5, 7, 23 };

    
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

    int32_t small_sieve_primes[25] =  { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    
    const uint32_t cache_bound = 100'000; //
    
    uint16_t i_p = 0;
    uint16_t i_L = 0;
    
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    
    mpz_cdiv_q( cmp_bound, bound, PL);
    
    std::vector< uint16_t > primes_lifting_L;
    
    uint64_t L_lift= 1;
    
    while( mpz_cmp_ui( cmp_bound, cache_bound ) > 0  )
    {
        if( small_sieve_primes[ i_p ] == L_distinct_primes[ i_L ] )
        {
            i_p++;
            if( i_L != L_len )
            {
                i_L++;
            }
        }
        else
        {
            L_lift *= small_sieve_primes[ i_p ];
            primes_lifting_L.push_back( small_sieve_primes[ i_p ] );
            std::cout << small_sieve_primes[ i_p ] << std::endl;
            mpz_mul_ui( PL, PL, small_sieve_primes[ i_p ] );
            mpz_cdiv_q( cmp_bound, bound, PL);
            i_p++;
        }
    }
    mpz_clear( cmp_bound );
    
    
    
    for( uint64_t m = 0; m < L_lift; m ++)
    {
        // for each r_star + m*L, we check that this number is not divisible by the small primes used to lift L
        bool enter_loop = true;
        for( auto p : primes_lifting_L )
        {
            enter_loop = ( enter_loop && ( mpz_divisible_ui_p( r_star, p ) == 0 ) );
        }
        
        // if r_star + m*L is not divisible by a lifted prime,
        // there are at most cache_bound modular exponentiations
        // to find potential CN
        if( enter_loop )
        {
            //set up the bit vector of sieving by the rest of the primes < append_bound here
            
            // this sets up n = P*(r_star + m*L)
            mpz_mul( n, P, r_star);
            
            
            // n is now created in an arithmetic progression of P*L*(lifted primes)
            while( mpz_cmp( n , bound ) < 0 )
            {
                // check bit vector before doing the modular exponentiation
                mpz_powm( base_2_fermat,  base2,  n, n); // 2^n mod n
                if( mpz_cmp( base_2_fermat, base2 ) == 0 )  // check if 2 = 2^n mod n
                {
                    mpz_powm( base_3_fermat,  base3,  n, n); // 3^n mod n
                    if( mpz_cmp( base_3_fermat, base3 ) == 0 )  // check if 3 = 3^n mod n
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
