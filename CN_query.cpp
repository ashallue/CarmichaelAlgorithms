//  We have work to do here....

// Answers the query "Is n a Carmichael number?" in 
// in deterministic polynomial time.
// See Algorithm 1 and Theorem 9 of draft
// Relies on GMP's B-PSW "primality" test
// A "no" answer most likely returns a Fermat b that witnessed in is composite
// A "yes" answer returns the complete prime factorization for verification with Korselt's
//        A "no" answer can be found with Korselt's criterion, too.  

// compiled with  g++ CN_search.cpp -lgmp -O3


#include <gmp.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>

    int main()
    {

    std::string string_number;
    std::cout << "This relies on GMP's B-PSW primality test and will answer the query:" << std::endl;
    std::cout << "Is n a Carmichael number?" << std::endl;
    std::cout << "Enter n: ";
    std::cin >> string_number; // Reads only the first word
        
	mpz_class n;
	n = string_number;

    if( mpz_probable_prime( n.get_mpz_t(), 10 ) > 0 )
    {
        std::cout << "GMP detects that n is a (probable) prime." << std::endl;
    }
    else
    {
        std::vector< mpz_class > current_composites;
        std::vector< mpz_class > next_composites;
        std::vector< mpz_class > primes;
        std::vector< mpz_class > fermat_base_powers;
        
        current_composite.push( n );
        
        mpz_class odd_part;
        odd_part = n - 1;
        // count powers of 2 dividing n-1:
        uint64_t pow_of_2 = (uint64_t) mpz_scan1( odd_part.get_mpz_t() , 0);
        // remove these powers of 2 to get the odd part:
        mpz_cdiv_q_2exp( odd_part.get_mpz_t(), n.get_mpz_t(), pow_of_2 );
        
        bool is_fermat_psp = true;
        
        mpz_class fermat_base = 1;
        
        while( !current_composites.empty() && is_fermat_psp )
        {
            
        }
    }
    

    
    

        
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    std::queue<uint64_t> R_composite_factors;
    std::vector<uint64_t> R_prime_factors;

    mpz_t fermat_result;
    mpz_t gcd_result;
    mpz_t r_factor;
    mpz_t base;
    mpz_t n;
    mpz_t odd_part;

    mpz_init( odd_part );
    mpz_init( n );
    mpz_init_set_ui( base, 2);
    mpz_init( r_factor );

    // SET P HERE
    // fermat_test is going to hold P until n is computed:
    //mpz_init_set_ui( fermat_result, 515410417841 );
    //uint64_t R = 1750661201;

    // below represents the least CN with 14 prime factors
    // where 6 of the primes are in R
    mpz_init_set_ui( fermat_result, 31792086053 );
    uint64_t R = 2757760839917;

    // set n = P*R, now fermat_result will be used in the obvious way
    mpz_mul_ui( n, fermat_result, R );

    // move R into mpz_t for prime test:
    mpz_init_set_ui( gcd_result, R );
    mpz_probab_prime_p( gcd_result, 0 ) == 0 ? R_composite_factors.push( R ) : R_prime_factors.push_back( R );

    // now using fermat_test to hold n-1
    mpz_sub_ui( fermat_result, n, 1);
    uint16_t pow_of_2 = (uint16_t) mpz_scan1( fermat_result, 0);
    mpz_fdiv_q_2exp( odd_part, fermat_result, pow_of_2 );

    uint64_t r64_factor;  			// holder for factors of R
    uint16_t start_size; 			// counter to empty queue
    uint16_t i = 0;					// counter for powers of 2
    bool is_fermat_psp = true;  	// initialized as true b/c n is a base-2 fermat psp

    // this while loops iterates over Fermat bases, starts with 2 and then odd integers
    mpz_powm( fermat_result, base, odd_part, n);
    mpz_set_ui( base, 1);
    while( !R_composite_factors.empty() && is_fermat_psp )
    {
        i = 0;	//counter for powers of 2
        // this while loop iterates over algebraic factors (b^(2^i*d) + 1)
        while( !R_composite_factors.empty() && i < pow_of_2 )
        {
            start_size = R_composite_factors.size();
            // this for loop iterates over composite factors
            for( int j = 0; j < start_size; j++ )
            {
                r64_factor = R_composite_factors.front();
                R_composite_factors.pop();
                mpz_set_ui( r_factor, r64_factor );
                mpz_add_ui( gcd_result, fermat_result, 1); // gcd_result hold b^(d*2^i) + 1
                mpz_gcd( gcd_result, gcd_result, r_factor);

                if( mpz_cmp(gcd_result, r_factor) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
                {
                    // this gcd split r_factor.  Letting g = gcd( gr, rf), then this splits rf into, g and rf/g.
                    r64_factor = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0 ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
                    mpz_divexact(gcd_result, r_factor, gcd_result );
                    r64_factor = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0  ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
                }
                else // r_factor was not factored, so put it back in the queue
                {
                    R_composite_factors.push( r64_factor );
                }
            }
            i++;
            mpz_powm_ui( fermat_result, fermat_result, 2, n );
        }
        
        std::cout << "We found " ;
        for( auto p : R_prime_factors )
        {
            std::cout << p << " " ;
        }
        std::cout << "as factors of R with the base " ;
        gmp_printf( "%Zd", base);
        std::cout << std::endl;
        
        mpz_add_ui( base, base, 2);
        mpz_powm( fermat_result, base, n, n);
        is_fermat_psp = ( mpz_cmp( fermat_result, base ) == 0 );
        mpz_powm( fermat_result, base, odd_part, n);
    }

    if( !is_fermat_psp )
    {
        std::cout << "The program stopped because a base " ;
        gmp_printf( "%Zd", base);
        std::cout << " Fermat test detected compositeness which cannot happen for CNs." << std::endl;
    }

    else if( R_composite_factors.empty( ) )
    {
        std::cout << "you now have a complete factorization - check with Korselt" << std::endl;
        std::cout << "the prime factors are " ;
        for( auto p : R_prime_factors )
        {
            std::cout << p << " " ;
        }
        std::cout << std::endl;
    }
        
    mpz_clear( fermat_result );
    mpz_clear( gcd_result );
    mpz_clear( r_factor );
    mpz_clear( base );
    mpz_clear( n );
    mpz_clear( odd_part );

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count() << "ms\n";

    }

