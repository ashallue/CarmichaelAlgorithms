// compiled with  g++ CN_search.cpp -lgmp -O3

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
    
    const uint32_t cache_bound = 100'000; //
    

    
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
        // for each r_star + m*L, we check that this number is not divisible by the small primes used to lift L
        bool enter_loop = true;
        for( auto p : primes_lifting_L )
        {
            enter_loop = ( enter_loop && ( mpz_divisible_ui_p( r_star, p ) == 0 ) );
        }
        
        // Check if we are in a valid spoke of a wheel
        // If so, sieve on that spoke
        if( enter_loop )
        {

            // set up the bit vector of sieving by the rest of the primes < append_bound here
            // the arithemtic progression below is (r_star + k* L_lift*L), where 0 <= k < cmp_bound
            // find k so that, r_star + k*(L_lift*L) = 0 mod p
            // so, we sieve out k = -r_star*(L_lift*L)^{-1} mod p
            spoke_sieve.reset();
            //  "20" has been hard-coded to this example.  Need to generalize be length of the primes array/vector
            mpz_t small_prime;
            mpz_init( small_prime );
            uint32_t k;
            uint16_t i_pp = i_p;
            while( i_pp < 20 )
            {
                // using "n" as a temp variable - it's set for its actual value immediately after this while loop ends
                mpz_set_ui( small_prime, small_sieve_primes[ i_pp ] );
                mpz_invert( n, PL, small_prime );  // n has PL^{-1} mod p
                mpz_neg( n, n);  // n has -(PL)^{-1} mod p
                mpz_mul( n, n, r_star ); // muliply by r_star
                mpz_mod( n, n, small_prime );  // reduce modulo p
                k = mpz_get_ui( n );
                
                // we just found the starting point for k, now sieve:
                while( k < cmp_bound32 )
                {
                    spoke_sieve[k] = 1;
                    k += small_sieve_primes[ i_pp ];
                }
                i_pp++;
            }
            mpz_clear( small_prime );

            
            // this sets up n = P*(r_star + m*L)
            mpz_mul( n, P, r_star);
            // n is now created in an arithmetic progression of P*L*(lifted primes)
            
           
            k = 0;
            while( mpz_cmp( n , bound ) < 0 )
            {
                if( spoke_sieve[k] == 0 )
                {
                    mpz_powm( base_2_fermat,  base2,  n, n); // 2^n mod n
                    if( mpz_cmp( base_2_fermat, base2 ) == 0 )  // check if 2 = 2^n mod n
                    {
                        mpz_powm( base_3_fermat,  base3,  n, n); // 3^n mod n
                        if( mpz_cmp( base_3_fermat, base3 ) == 0 )  // check if 3 = 3^n mod n
                        {
                            // n is now base 2 and a base 3 Fermat psp
                            // invoke CN factorization algorithm here

                            
                            std::cout << "The value for m is " << m << std::endl;
                            std::cout << "The value for k is " << k << std::endl;
                            std::cout << "The value of PL is " ;
                            gmp_printf( "%Zd", PL);
                            std::cout << std::endl;
                            std::cout << "The value of r_star is " ;
                            gmp_printf( "%Zd", r_star);
                            std::cout << std::endl;

                            std::cout << "n = " ;
                            gmp_printf( "%Zd", n);
                            std::cout << " is a base-2 and base-3 Fermat psp." << std::endl;
                            
                            
                        }
                    }
                }
                k++;
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
