#include "Preproduct.h"
#include "rollsieve.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <queue>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>
#include <cstddef>
#include <boost/dynamic_bitset.hpp>

static_assert(sizeof(unsigned long) == 8, "unsigned long must be 8 bytes.  needed for mpz's unsigned longs to take 64 bit inputs in various calls.  LP64 model needed ");

// preprocessing flag.  If enabled, bounds and some other assumptions will change for testing purposes
#define TEST

Preproduct::Preproduct()
{
    mpz_init( P );
    mpz_init( L );
    mpz_init_set_ui( BOUND, 10 );

    #ifdef TEST
        // bound for testing
        mpz_pow_ui( BOUND, BOUND, 15 );
    #else
        // bound for full computation
        mpz_pow_ui( BOUND, BOUND, 24 );
    #endif
    
}

Preproduct::~Preproduct()
{
    mpz_clear( P );
    mpz_clear( L );
    mpz_clear( BOUND );
}

// assumes valid inputs: 
// 1) does not check that init_preproduct is cylic 
// 2) does not check that init_LofP is actually CarmichaelLambda( init_preproduct )
// intended use is initializing from precomputation which only generates valid inputs
// factors P using trial division because the intended use case is relatively small initializing preproducts
// could consider another version that passes array/vector with primes
// could also consider faster factorization algorithms
void Preproduct::initializing( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound )
{
    uint64_t temp;
        
    //set info for P:
    mpz_set_ui( P, init_preproduct );
    
    // clear vector and fill
    P_primes.clear();
    if( init_preproduct != 1 )
    {
        temp = 1;
        while( temp*temp < init_preproduct )
        {
            temp+= 2;
            if( ( init_preproduct % temp ) == 0 )
            {
                init_preproduct = init_preproduct / temp;
                P_primes.push_back(temp);
            }
        }
        // possible remaining factor
        if( init_preproduct > 1 )
        {
            P_primes.push_back( init_preproduct );
        }
    }
    
    //set info for L:
    mpz_set_ui( L, init_LofP );
    
    // set append_bound
    append_bound = init_append_bound;

}


// assumes prime is admissible to PP
void Preproduct::appending(Preproduct& PP, uint64_t prime )
{
    // compute new P value:
    mpz_set( P, PP.P );
    mpz_mul_ui( P, P, prime );
    
    // set the vector
    P_primes.clear();
    P_primes = PP.P_primes;
    P_primes.push_back( prime );

    // compute L
    mpz_set( L, PP.L );
    mpz_lcm_ui( L, L, prime - 1 );

    append_bound = prime;
}


bool Preproduct::is_admissible_modchecks( uint64_t prime_to_append )
{
    bool return_val = true;
    
    // can be done with std::all_of or std::any_of or std::none_of
    uint16_t i = 0;
    while( return_val && i < P_primes.size() )
    {
        return_val = return_val && (1 != ( prime_to_append % P_primes[i] ) );
        i++;
    }
 
    return return_val;
}

void Preproduct::complete_tabulation( std::string cars_file )
{
    #ifdef TEST
        // looking for cars that are duplicates
        if(mpz_cmp_ui(P, 1) == 0)
        {
            std::cout << "Calling complete_tabulation with P = 1\n ";
        }
    #endif
    
    
    // Unanswered question that we need to answer before production:
    // Do we need to consider if P is, itself, a CN in this?  If so, where is that done?  right here?
    // Answer: yes. Call appending_is_CN with an empty vector. Check that P.p_primes.size() > 2
    if( P_primes.size() > 2 ) 
    {   std::vector< uint64_t> empty;
        appending_is_CN( empty, cars_file);
    }
    
    // rule to be set later
    // We show that if PL^2 > B then rule should be true
    // however, assumes the cost of the else block is porportional to the cost of generating the primes
    // this is almost certain too small of an estimate.
    // something like PL^3 > B might be justified
    // some analytic work required but not needed for correctness
    
    mpz_t to_recurse_or_not_to_recurse;
    mpz_init_set(  to_recurse_or_not_to_recurse, P);
    mpz_mul(  to_recurse_or_not_to_recurse,  to_recurse_or_not_to_recurse, L);
    mpz_mul(  to_recurse_or_not_to_recurse,  to_recurse_or_not_to_recurse, L);
    bool rule = ( mpz_cmp( to_recurse_or_not_to_recurse , BOUND ) >= 0 );
   
    if( rule )
    {
        #ifdef TEST
            if(mpz_cmp_ui(P, 79003) == 0) std::cout << "P = 79003, calling CN_search\n";
        #endif
        
        CN_search( cars_file );
    }
    else
    {
        #ifdef TEST
            if(mpz_cmp_ui(P, 79003) == 0) std::cout << "P = 79003, recursing\n";
        #endif
        
        
        // in the below, let p be the largest prime dividing P
        // by taking this route, we need to do up to three things:
        // 1 - find the single prime q so that Pq is a CN
        // 2 - for q in ( append_bound, (B/P)^(1/3) ), Pq is small enough to recurse on the preproduct Pq
        // 3 - for q in  ( (B/P)^(1/3) ,  (B/P)^(1/2) ), Pq can only have a single prime appended

        // tabulation assumes all CN of the form P*q*r with P < X have already been found
        // This creates 2 special cases for the above cases:
        //      A - if P < X, then only need case 2 - (case 1 and case 3 involve appending exactly 1 or exactly 2 primes)
        //      B - if P/p < X, then we only need case 2 and case 3 (if a single prime is found then (P/p)*p*q is a CN and P/p < X , so it is duplicated)

        // X set to 120m for production run, because of existing tabulation of the small case
        // X set to 1 for testing, because testing focused on the large case
        #ifndef TEST
            const uint64_t X = 120'000'000;
        #else
            const uint64_t X = 1;
        #endif
        
        mpz_t BoverP;
        mpz_t case_bound;
        mpz_init( case_bound );
        mpz_init( BoverP);
        mpz_cdiv_q( BoverP, BOUND, P );
                
        // this is case 1
        // this is only invoked if P > X * P_primes.back()
        // the count of prime factors in P needs to be greater than 1
        // if X > B^(1/3), then the count needs to be greater than 2
        // because all CN with exactly three prime factors has already been found
        mpz_set_ui( case_bound, X);

        // only multiply if there is a prime to multiply by
        if( P_primes.size() > 0 )
        {
            mpz_mul_ui( case_bound, case_bound, P_primes.back() );
        }

        // for testing complete with one prime if preproduct has 2 or more primes
        // for production, the 3-carmichael case entirely small, so large case restricted to 4 prime factors or more
        #ifdef TEST
            if( P_primes.size() > 1 && mpz_cmp( P, case_bound ) > 0 )
        #else
            if( P_primes.size() > 2 && mpz_cmp( P, case_bound ) > 0 )
        #endif
        {
            completing_with_exactly_one_prime( cars_file );
        }
        
        // this is the start of cases 2 and 3: they share the incremental sieve and form a preproduct Pq
        Preproduct Pq;
        Rollsieve r( append_bound + 1 );
       
        // this is the start of case 2.  case2_bound stores (B/P)^{1/3}
        mpz_root( case_bound, BoverP, 3);
        uint64_t case2_bound = mpz_get_ui( case_bound );
        
        
        uint64_t q = r.nextprime();
        while( q < case2_bound )
        {
            
            if( is_admissible_modchecks( q ) )
            {
                
                Pq.appending( *this, q ) ;
                Pq.complete_tabulation( cars_file );
            }
            q = r.nextprime();
        }
        
        // this is the start of caes 3
        // first we check that it needs to be done
        // see comments on case 1 regarding P_primes.size()
        // since this will be invoked on a preproduct Pq
        // the inequality on prime counts lowered by 1
        #ifdef TEST
            if( P_primes.size() > 0 && mpz_cmp_ui( P, X ) > 0 )
        #else
                
            if( P_primes.size() > 1 && mpz_cmp_ui( P, X ) > 0 )
        #endif
        {
            
            mpz_sqrt( case_bound, BoverP );
            uint64_t case3_bound = mpz_get_ui( case_bound );
            
            while( q < case3_bound )
            {
                if( is_admissible_modchecks( q ) )
                {         
                    #ifdef TEST
                    if(mpz_cmp_ui(P, 79003) == 0 && q < 5000) std::cout << "P = 79003, case 3, q = " << q << "\n";
                    #endif

                    // append q to the current object
                    Pq.appending( *this, q ) ;

                    // if Pq itself is Carmichael, write it to file
                    std::vector< uint64_t> empty;
                    appending_is_CN( empty, cars_file);

                    // then attempt to complete with one more prime past Pq
                    Pq.completing_with_exactly_one_prime( cars_file );
                }
                q = r.nextprime();
            }
        }
        
        mpz_clear( BoverP );
        mpz_clear( case_bound );
    }
    
    mpz_clear( to_recurse_or_not_to_recurse );
}



// Do not call this method on a preproduct of the form (1, 1, b)
// Calling this method on (1, 1, b) results in a lienar search up to B
// and *will* return prime numbers as Carmichael numbers because 1*p passes the Korselt check
// the rule used in complete tabulation should prevent (1,1,b) from being used to call this
void Preproduct::CN_search( std::string cars_file )
{
    mpz_t r_star;
    mpz_init( r_star );
    mpz_invert(r_star, P, L);
    
    mpz_t base2;
    mpz_t base3;
    mpz_init_set_ui( base2, 2 );
    mpz_init_set_ui( base3, 3 );
    
    mpz_t fermat_result;
    mpz_init( fermat_result );
    
    mpz_t n;
    mpz_init(n);

    // n = P*( r_star + k*L ) = P*r_star * k*(P*L), an arithmetic progression with common difference PL
    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );
    
    mpz_t small_prime;
    mpz_init( small_prime );
    
    // position i has the truth value of the statement "(2i + 3) is prime"
    const uint32_t bitset_size = 512;
    std::bitset<bitset_size> small_primes{"00110010100000100100010010010100000010010010100010000100010100000000010110100000010110100000010000110110000110000010000100000010100010100100010100100100010000100010000100010010100000110010010110000100000110100100110010010000100110010010000100100000000110000010010100010100010000010110100010010100110000110000100010100010010100100100010010110000100100000010110100000010000110100110010010010000110010110100000100000110110000110010010100100110000110010100000010110110100010010100110100110010010110100110010110110111"};
    
    mpz_t R;
    mpz_init( R );                              // using R as a temp variable here
    mpz_cdiv_q( R, BOUND, PL);
    
    uint64_t cmp_bound64 = mpz_get_ui( R );     // done using R as a temp variable
    boost::dynamic_bitset<> spoke_sieve( cmp_bound64 );
    spoke_sieve.reset();
    
    // sieve
    uint16_t prime_index = 0;
    while( prime_index < bitset_size )
    {
        if( small_primes[ prime_index ] )
        {
            uint32_t p = 2*prime_index + 3;
            // we do not want R to be divisible by primes less than append_bound
            // we can also use any prime inadmissible to P
            if( p <= append_bound || !is_admissible_modchecks( p ) )
            {
                mpz_set_ui( small_prime, p );
                if( mpz_invert( n, L, small_prime ) )   // n has L^{-1} mod p, if the inverse exists
                {
                    mpz_neg( n, n);                     // n has -(L)^{-1} mod p
                    mpz_mul( n, n, r_star );            // muliply by r_star
                    mpz_mod( n, n, small_prime );       // reduce modulo p
                    uint64_t k = mpz_get_ui( n );       // k has the starting point for sieve
                    while( k < cmp_bound64 )
                    {
                        spoke_sieve[k] = 1;
                        k += p;
                    }
                }
            }
        }
        prime_index++;
    }
    // sieving now done
    
    mpz_mul( n, P, r_star);     // n is now iniitalized as P * r_star + 0*PL
    uint32_t k = 0;             // and so k has the value 0
    while( mpz_cmp( n , BOUND ) < 0 )
    {
        if( spoke_sieve[k] == 0 )
        {
            mpz_powm( fermat_result,  base2,  n, n);        // 2^n mod n
            if( mpz_cmp( fermat_result, base2 ) == 0 )      // check if 2 = 2^n mod n
            {
                mpz_powm( fermat_result,  base3,  n, n);    // 3^n mod n
                if( mpz_cmp( fermat_result, base3 ) == 0 )  // check if 3 = 3^n mod n
                {
                    mpz_divexact( R, n, P);
                    CN_factorization( n, R, cars_file  );
                    // gmp_printf( "n = %Zd = %Zd * %Zd is a base-2 and base-3 Fermat psp. \n", n, P, R);
                }
            }
        }
        k++;
        mpz_add( n, n, PL);
    }
    
    mpz_clear( R );
    mpz_clear( small_prime );
    mpz_clear( r_star );
    mpz_clear( n );
    mpz_clear( PL );
    mpz_clear( base2 );
    mpz_clear( base3 );
    mpz_clear( fermat_result );

}


/* Depends on primes_to_append being a vector of true primes.
    
*/
bool Preproduct::appending_is_CN( std::vector< uint64_t >&  primes_to_append, std::string cars_file )
{
    // file object for storing the carmichael numbers
    FILE* cars_output;
    //char filename[100];
    //cars_file.copy(filename, cars_file.length());
    const char* filename;
    filename = cars_file.c_str();
    
    cars_output = fopen (filename,"a");
    
    mpz_t P_temp;
    mpz_t L_temp;

    mpz_init_set( P_temp, P );
    mpz_init_set( L_temp, L );

    // "update" P and L for each prime
    for( auto app_prime : primes_to_append )
    {
        // update the preproduct by this prime
        mpz_mul_ui( P_temp, P_temp, app_prime );
        uint64_t temp = app_prime - 1;
        mpz_lcm_ui( L_temp, L_temp, temp );
    }
    
    mpz_sub_ui( P_temp, P_temp, 1);
    bool return_val =  mpz_divisible_p( P_temp, L_temp );
    mpz_add_ui( P_temp, P_temp, 1);

    // if return_val, then L | n-1, making n Carmichael.  Write n to file
    // we also check a boundedness condition to ensure n is within upper bound
    if( return_val && mpz_cmp( P_temp, BOUND ) < 0 )
    {
        #ifdef TEST
            // looking for cars that are duplicates
            if(mpz_cmp_ui(P_temp, 41041) == 0)
            {
                std::cout << "41041 found with preproduct: ";
                gmp_printf("%Zd\n", P);
            }
        #endif
        
        gmp_fprintf(cars_output, "%Zd\n", P_temp );
    }

    mpz_clear( P_temp );
    mpz_clear( L_temp );
    
    // close file
    fclose (cars_output);

    return return_val;
}

// Check that lambda(P) divides (P-1)
bool Preproduct::is_CN( )
{
    bool return_val;
    mpz_sub_ui( P, P, 1);
    return_val = mpz_divisible_p( P, L );
    mpz_add_ui( P, P, 1);
    return return_val;
}


// see section 5.3 of ANTS 2024 work and the corresponding implementation:
// https://github.com/ashallue/tabulate_car/blob/master/LargePreproduct.cpp#L439C1-L500C2
void Preproduct::completing_with_exactly_one_prime( std::string cars_file )
{
    #ifdef TEST
        // looking for cars that are duplicates
        if(mpz_cmp_ui(P, 115921) == 0)
        {
            std::cout << "P = 115921 in completing_with_exactly_one_prime\n";
        }
        if(mpz_cmp_ui(P, 481) == 0)
        {
            std::cout << "P = 481, completing with exactly one prime\n";
        }
    #endif
    
    std::vector <uint64_t> the_prime_factor;
    uint64_t prime_factor;
    
    // construct r* = P^{-1} mod L, the least positive residue
    mpz_t r_star;
    mpz_init( r_star );
    mpz_invert( r_star, P, L);
    
    // construct g = gcd( r* - 1, L),  the scaling factor
    mpz_t g;
    mpz_init( g );
    mpz_sub_ui( g, r_star, 1);  // it holds r^{\star} - 1
    mpz_gcd( g, g, L );    // g = gcd ( r^{\star} - 1, L )
    
    // construct script_P = (P - 1)/g.  
    // Note that r-1 = r*-1 mod L and r-1 | P-1 by assumption, so P-1 is div by g
    mpz_t script_P;
    mpz_init_set( script_P, P );
    mpz_sub_ui( script_P, script_P, 1 );
    mpz_divexact( script_P, script_P, g);
        
    // set div_bound1.  This is sec 5.3.1 of the Advances paper
    // script_R + k*script_L < sqrt( script_P )  iff
    // g*(script_R + k*script_L) + 1 < g*sqrt( script_P ) + 1  iff
    // r* + k*L < g*sqrt( script_P ) + 1
    mpz_t div_bound1;
    mpz_init_set( div_bound1, script_P );
    mpz_mul( div_bound1, div_bound1, g );
    mpz_mul( div_bound1, div_bound1, g );
    mpz_sqrt( div_bound1, div_bound1 );
    // [Andrew]     Does this round down in a way we don't want? sqrt does truncate
    // [Jonathan]   It is fine to truncate.
    //              We use less than or equal to account for script_P being a perfect square
    //              otherwise this is not a big deal
    mpz_add_ui( div_bound1, div_bound1, 1);
    
    // set div_bound2
    // P( r* + k*L ) <  B
    // r* + k*L < B/P
    mpz_t div_bound2;
    mpz_init_set( div_bound2, BOUND );
    // [Andrew] again, rounding down with strict inequality might be bad.  Case 1 has non-strict inequality at least.
    // [Jonathan]  cdiv is the ceiling division.  This does not round down.
    mpz_cdiv_q( div_bound2, div_bound2, P);

    // returns true if div_bound1 > div_bound2.
    // This corresponds to the beginning of section 5.2 from Advances
    // check arithmetic progression r* + kL up to B/P
    if( mpz_cmp( div_bound1, div_bound2 ) > 0 )
    {
        // div_bound2 = B/P is the relevant bound
        // we do not need to look for the large divisor (see 5.3.2 Case 2)
        // this path is essentially CN_search but with no sieving and no CN_factorization needed
        while( mpz_cmp( r_star, div_bound2) <= 0 )
        {
            // primality check on r* + kL
            if( mpz_probab_prime_p( r_star, 0) != 0 )
            {
                prime_factor = mpz_get_ui( r_star );
                if( prime_factor > append_bound )
                {
                    the_prime_factor.clear();
                    the_prime_factor.push_back( prime_factor );

                    #ifdef TEST
                        // looking for cars that are duplicates
                        if(mpz_cmp_ui(P, 481) == 0 && prime_factor == 241)
                        {
                            std::cout << "P = 481, case 1\n";
                        }
                    #endif
                    // call helper function to check if Pr is Carmichael.  That function will write to file
                    appending_is_CN( the_prime_factor , cars_file );
                }
            }
            // next iteration: add L
            mpz_add( r_star, r_star, L );
        }
    }
    // In this case div_bound1 is smaller, so we search for r* + kL < g*sqrt (script_P) + 1
    else
    {
        // we need script_L
        mpz_t script_L;
        mpz_init( script_L );
        mpz_divexact( script_L, L, g );

        // we need to use r_star again [Andrew] - I don't get this
        // prev code was r2 = ( inv128( r1, L1) * scriptP ) % L1 where r1 = (Pqinv - 1)/g
        // Here this becomes r1 = (r* - 1)/g, r2 = (r1^{-1} mod scriptL) * scriptP mod scriptL
        mpz_t r2;
        mpz_init( r2 );

        // r2 is going to hold r1 = (r* - 1)/g
        mpz_set( r2, r_star );
        mpz_sub_ui( r2, r2, 1);
        mpz_divexact( r2, r2, g);

        // r1 computed, now set r2 = script_P*(r1)^{-1} mod script_L
        mpz_invert( r2, r2, script_L );
        mpz_mul( r2, r2, script_P );
        mpz_mod( r2, r2, script_L );

        // Following section 5.3.1 of Advances, scriptR + k scriptL equiv to arithmetic progression r* + kL
        // div_bound1 reflects the transformation: multiply by g and add 1.
        while( mpz_cmp( r_star, div_bound1) <= 0 )
        {
            // primality check.
            if( mpz_probab_prime_p( r_star, 0) > 0 )
            {
                prime_factor = mpz_get_ui( r_star );
                if( prime_factor > append_bound )
                {
                    the_prime_factor.clear();
                    the_prime_factor.push_back( prime_factor );

                    #ifdef TEST
                        // looking for cars that are duplicates
                        if(mpz_cmp_ui(P, 481) == 0 && prime_factor == 241)
                        {
                            std::cout << "P = 481, case 2\n";
                        }
                    #endif
                     // call helper function to check if Pr is Carmichael.  That function will write to file
                    appending_is_CN( the_prime_factor , cars_file );
                }
            }
            // next iteration
            mpz_add( r_star, r_star, L );
        }
        
        // can do a computation to find a possilbe k > 0 so and set
        // r2 = r2 + k*script_L
        // see this done in the previous version:
        // https://github.com/ashallue/tabulate_car/blob/master/LargePreproduct.cpp#L479-L484

        // resetting div_bound1 to now hold floor sqrt script_P
        mpz_sqrt( div_bound1, script_P );

        #ifdef TEST
        if (mpz_cmp_ui(P, 481) == 0)
        {
            gmp_printf("div_bound1 = %Zd, script_P = %Zd\n", div_bound1, script_P);
        }
        #endif
        // check arith progression r2 + k*scriptL < sqrt(script_P)
        while( mpz_cmp( r2, div_bound1) < 0 )
        {
            #ifdef TEST
            if (mpz_cmp_ui(P, 481) == 0)
            {
                std::cout << "r2 = ";
                gmp_printf("%Zd\n", r2 );
            }
            #endif
            if( mpz_divisible_p( script_P, r2 ) )
            {
                // we are done with r_star, using it as storage for R.  Possible prime is r2*g + 1
                mpz_divexact( r_star, script_P, r2);
                mpz_mul( r_star, r_star, g);
                mpz_add_ui( r_star, r_star, 1 );
                if( mpz_probab_prime_p( r_star,  0 ) > 0 )
                {
                    prime_factor = mpz_get_ui( r_star );
                    if( prime_factor > append_bound )
                    {
                        the_prime_factor.clear();
                        the_prime_factor.push_back( prime_factor );

                        #ifdef TEST
                            // looking for cars that are duplicates
                            if(mpz_cmp_ui(P, 481) == 0 && prime_factor == 241)
                            {
                                std::cout << "P = 481, case 3\n";
                            }
                        #endif
                        
                        appending_is_CN( the_prime_factor , cars_file );
                    }
                }
            }
            // add script_L to continue the progression
            mpz_add( r2, r2, script_L );
        }

        mpz_clear( r2 );
        mpz_clear( script_L );
    }
    
    mpz_clear( div_bound2 );
    mpz_clear( div_bound1 );
    mpz_clear( script_P );
    mpz_clear( g );
    mpz_clear( r_star );
}


// R is passed as mpz_t type but we assume is that R < 2^64 because it is immediately written put into a uint64_t type
// minor improvements:
//  1) break earlier with append bound
//  2) b^n mod n is, in essence, computed twice for each b - only compute it once
// Return type is bool.  False means n failed a fermat test, so is not Carmichael.
bool Preproduct::CN_factorization( mpz_t& n, mpz_t& R, std::string cars_file )
{
    std::queue<uint64_t> R_composite_factors;
    std::vector<uint64_t> R_prime_factors;
    
    uint64_t r64_factor =  mpz_get_ui( R );
    mpz_probab_prime_p( R, 0 ) == 0 ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
    
    mpz_t fermat_result;
    mpz_t odd_part;
    mpz_t base;
    mpz_t gcd_result;
    
    mpz_init( gcd_result );
    mpz_init_set_ui( base, 2);
    
    mpz_init( odd_part );
    mpz_init_set( fermat_result, n );
    // now using fermat_test to hold n-1
    mpz_sub_ui( fermat_result, fermat_result, 1 );
    uint16_t pow_of_2 = (uint16_t) mpz_scan1( fermat_result, 0 );
    mpz_fdiv_q_2exp( odd_part, fermat_result, pow_of_2 );
    
    uint16_t start_size;            // counter to empty queue
    uint16_t i = 0;                 // counter for powers of 2
    bool is_fermat_psp = true;      // initialized as true b/c n is a base-2 fermat psp

    // this while loops iterates over Fermat bases, starts with 2 and then odd integers
    mpz_powm( fermat_result, base, odd_part, n);
    mpz_set_ui( base, 1);
    while( !R_composite_factors.empty() && is_fermat_psp )
    {
        i = 0;    //counter for powers of 2
        // this while loop iterates over algebraic factors (b^(d * 2^i) + 1)
        // early abort?  no need to check all factors
        while( !R_composite_factors.empty() && i < pow_of_2 )
        {
            start_size = R_composite_factors.size();
            // this for loop iterates over composite factors
            for( int j = 0; j < start_size; j++ )
            {
                r64_factor = R_composite_factors.front();
                R_composite_factors.pop();
                mpz_set_ui( R, r64_factor );
                mpz_add_ui( gcd_result, fermat_result, 1); // gcd_result hold b^(d*2^i) + 1
                mpz_gcd( gcd_result, gcd_result, R);

                if( mpz_cmp(gcd_result, R) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
                {
                    // this gcd split r_factor.  Letting g = gcd( gr, rf), then this splits rf into, g and rf/g.
                    // could break here if rf or rf/g < append_bound (exiting early if not)
                    r64_factor = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0 ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
                    mpz_divexact(gcd_result, R, gcd_result );
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

        // get new base
        mpz_add_ui( base, base, 2);
        mpz_powm( fermat_result, base, n, n);
        is_fermat_psp = ( mpz_cmp( fermat_result, base ) == 0 );
        mpz_powm( fermat_result, base, odd_part, n);
    }

    if( R_composite_factors.empty( ) && is_fermat_psp )
    {
        // could grab the minimum of R_prime_factors but we are sorting for the output to be correct
        std::sort ( R_prime_factors.begin(), R_prime_factors.end() );

        if( R_prime_factors[0] > append_bound )
        {
            appending_is_CN( R_prime_factors , cars_file);
        }
    }
        
    mpz_clear( fermat_result );
    mpz_clear( gcd_result );
    mpz_clear( base );
    mpz_clear( odd_part );
    
    return is_fermat_psp;
}





