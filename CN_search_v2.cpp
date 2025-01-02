// compiled with  g++ CN_search.cpp -lgmp -O3

// wheel
// no spoke sieving

#include <gmp.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>
#include <boost/dynamic_bitset.hpp>

int main()
{
    // there are many hard-coded things in the below that do not make this easy to change.
    // in particular the "small primes" vector is missing primes that divide L
    // here is an attempt at tracking all the places something is hard-coded
    // The bound - see lines 34-36
    // P - see line 40
    // CarmichaelLambda(P) = L - See line 44 - it is not check that this quantity is, indeed, L(P)
    // the primes dividing L - See line 45 and 46
    // the small primes that aren't in L - see line 78 and 135 were a hard-coded "20" also appears
    // the cache_bound - a rough guess to make sure the dynamic bitset fits in L1 cache
    
    
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

    // we will need to remove from this vector/array the primes dividing L
    // this is hard-coded for this example right now
    int32_t small_sieve_primes[20] =  { 11, 13, 17, 19, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    
    const uint32_t cache_bound = 200'000; //
   
    
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    
    mpz_cdiv_q( cmp_bound, bound, PL);
    
    std::vector< uint16_t > primes_lifting_L;
    
    uint64_t L_lift= 1;
    
    uint16_t i_p = 0;
    
    while( mpz_cmp_ui( cmp_bound, cache_bound ) > 0  )
    {
        L_lift *= small_sieve_primes[ i_p ];
        primes_lifting_L.push_back( small_sieve_primes[ i_p ] );
        mpz_mul_ui( PL, PL, small_sieve_primes[ i_p ] );
        mpz_cdiv_q( cmp_bound, bound, PL);
        i_p++;
    }

    // cmp_bound needs to be put in a uint32_t
    uint32_t cmp_bound32 = mpz_get_ui( cmp_bound );
    boost::dynamic_bitset<> spoke_sieve( cmp_bound32 );
    
    std::cout << "The number of residues per spoke is " << cmp_bound32 << std::endl;
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
        
        if( enter_loop )
        {
            mpz_mul( n, P, r_star);
            while( mpz_cmp( n , bound ) < 0 )
            {
                mpz_powm( base_2_fermat,  base2,  n, n); // 2^n mod n
                if( mpz_cmp( base_2_fermat, base2 ) == 0 )  // check if 2 = 2^n mod n
                {
                    mpz_powm( base_3_fermat,  base3,  n, n); // 3^n mod n
                    if( mpz_cmp( base_3_fermat, base3 ) == 0 )  // check if 3 = 3^n mod n
                    {
                        // CN check here.
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
       
    mpz_clear( cmp_bound );
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
