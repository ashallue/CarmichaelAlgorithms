// compiled with  g++ CN_search.cpp -lgmp -O3

#include <gmp.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>

int main()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    // This program assumes n = P*R satisifies
    // 2^n = 2 mod n and 3^n = 3 mod n
    
    mpz_t P;
    mpz_init_set_ui( P, 515410417841 );

    uint64_t R = 197701945121;
    mpz_t gmpR;
    mpz_init_set_ui( gmpR, R);
    
    /*
     all work with this P value
     1750661201
     7166791361
     20531787761
     31883253761
     197701945121
     292890464561
     305673542561
     501352298561
     654409588961
     654775316561
     757543033361
     1085677545521
     1303915724561
     1417671498161
     */

    
    mpz_t n;
    mpz_init_set( n, P);
    mpz_mul_ui( n, P, R );
    
    mpz_t nm1;
    mpz_init( nm1 );
    mpz_sub_ui( nm1, n, 1);
    uint16_t pow_of_2 = (uint16_t) mpz_scan1( nm1, 0);
    
    mpz_t odd_part;
    mpz_init( odd_part );
    mpz_fdiv_q_2exp( odd_part, nm1, pow_of_2 );
    
    mpz_t base;
    mpz_init_set_ui( base, 2);
    
    mpz_t fermat_result;
    mpz_init( fermat_result );
    mpz_t fermat_result_p1;
    mpz_init( fermat_result_p1 );
    
   
    mpz_t gcd_result;
    mpz_init( gcd_result );
    mpz_t r_factor;
    mpz_init( r_factor );
    
    std::cout << "n = " ;
    gmp_printf( "%Zd", n);
    std::cout << " is a base-2 and base-3 Fermat psp." << std::endl;
    std::cout << "n - 1 = " ;
    gmp_printf( "%Zd", nm1);
    std::cout << "." << std::endl;
    std::cout << "P = " ;
    gmp_printf( "%Zd", P);
    std::cout << "." << std::endl;
    std::cout << "Odd part of n-1 is " ;
    gmp_printf( "%Zd", odd_part);
    std::cout << "." << std::endl;
    std::cout << "n-1 has " << pow_of_2 << " powers of 2"<< std::endl;
    std::cout << "R = " ;
    gmp_printf( "%Zd", gmpR);
    std::cout << " = " << R <<  std::endl;

    std::queue<uint64_t> R_composite_factors;
    std::vector<uint64_t> R_prime_factors;
    
    bool is_fermat_psp = true;
    
    if( mpz_probab_prime_p( gmpR, 0 ) == 0 )
    {
        std::cout << "R is not prime" << std::endl;
        R_composite_factors.push( R );
    }
    else
    {
        std::cout << "R is prime." << std::endl;
    }

    uint64_t temp;
       
    uint16_t start_size;
    mpz_powm( fermat_result, base, odd_part, n);
    uint16_t i = 0;
    while( !R_composite_factors.empty() && i < pow_of_2 )
    {
        start_size = R_composite_factors.size();
        for( int j = 0; j < start_size; j++ )
        {
            temp = R_composite_factors.front();
            R_composite_factors.pop();
            mpz_set_ui( r_factor, temp );
            
            mpz_add_ui( fermat_result_p1, fermat_result, 1);
            mpz_gcd( gcd_result, fermat_result_p1, r_factor);
            
            if( mpz_cmp(gcd_result, r_factor) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
            {
                uint64_t temp = mpz_get_ui( gcd_result );
                mpz_probab_prime_p( gcd_result, 0 ) == 0 ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
                mpz_divexact(gcd_result, r_factor, gcd_result );
                temp = mpz_get_ui( gcd_result );
                mpz_probab_prime_p( gcd_result, 0 ) == 0  ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
            }
            else // r_factor was not factored, so it is prime or composite
            {
                R_composite_factors.push( temp );
            }
        }
        i++;
        mpz_powm_ui( fermat_result, fermat_result, 2, n );
        std::cout << R_composite_factors.size() << " " << R_prime_factors.size() << std::endl;
    }
    
    std::cout << "We found " ;
    for( auto p : R_prime_factors )
    {
        std::cout << p << " " ;
    }
    std::cout << "as factors of R with the base 2" << std::endl;
    
    
    
    
    mpz_set_ui( base, 3);
    mpz_powm( fermat_result, base, n, n);
    is_fermat_psp = ( mpz_cmp( fermat_result, base ) == 0 );
        
    while( !R_composite_factors.empty() && is_fermat_psp )
    {
        mpz_powm( fermat_result, base, odd_part, n);
        i = 0;
        while( !R_composite_factors.empty() && i < pow_of_2 )
        {
            start_size = R_composite_factors.size();
            for( int j = 0; j < start_size; j++ )
            {
                temp = R_composite_factors.front();
                R_composite_factors.pop();
                mpz_set_ui( r_factor, temp );
                
                mpz_add_ui( fermat_result_p1, fermat_result, 1);
                mpz_gcd( gcd_result, fermat_result_p1, r_factor);
                
                if( mpz_cmp(gcd_result, r_factor) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
                {
                    uint64_t temp = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0 ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
                    mpz_divexact(gcd_result, r_factor, gcd_result );
                    temp = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0  ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
                }
                else // r_factor was not factored, so put it back in the queue
                {
                    R_composite_factors.push( temp );
                }
            }
            i++;
            mpz_powm_ui( fermat_result, fermat_result, 2, n );
            std::cout << R_composite_factors.size() << " " << R_prime_factors.size() << std::endl;
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
    }
    
    if( !is_fermat_psp )
    {
        std::cout << "The program stopped because a simple base " ;
        gmp_printf( "%Zd", base);
        std::cout << " Fermat test detected compositeness (which cannot happen for CNs)." << std::endl;
    }
    
    if( R_composite_factors.empty( ) )
    {
        std::cout << "you now have a complete factorization - check with Korselt" << std::endl;
    }
    
    
        
    mpz_clear( r_factor );
    mpz_clear( gcd_result );
    mpz_clear( gmpR ) ;
    mpz_clear( odd_part );
    mpz_clear( nm1 );
    mpz_clear( n );
    mpz_clear( P );
    
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count() << "ms\n";



}

