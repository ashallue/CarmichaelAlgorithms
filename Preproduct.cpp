#include "Preproduct.h"
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

// redo these if necessary

// hard code sqrt( B )
// we have choosen B = 10^24
#define SQRT_BOUND 1'000'000'000'000
// largest prime <= sqrt( B / X ) = 10^8
// because X = 10^8
#define DEFAULT_MAX_PRIME_BOUND 100'000'000
// there are 5761455 primes less than 10^8
#define PRIME_COUNT 5761455

Preproduct::Preproduct()
{
    mpz_init( P ) ;
    mpz_init( L ) ;
    mpz_init_set_ui( BOUND, 10 );
    mpz_pow_ui( BOUND, BOUND, 24 );
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
// uses trial division because the intended use case is relatively small initializing preproducts
// could consider another version that passes in arrays of factors
// precomputation.cpp would need to be re-written
void Preproduct::initializing( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound )
{
    mpz_set_ui( P, init_preproduct );
    mpz_set_ui( L, init_LofP );
    append_bound = init_append_bound;
    P_len = 0;

       
    uint64_t temp;
    // set primes array for P
    if( init_preproduct != 1 )
    {
        temp = 1;
        while( temp*temp < init_preproduct )
        {
            temp+= 2;
            if( ( init_preproduct % temp ) == 0 )
            {
                init_preproduct = init_preproduct / temp;
                P_primes[ P_len ] = temp;
                P_len++;
            }
        }
        // possible remaining factor
        if( init_preproduct > 1 )
        {
            P_primes[ P_len ] = init_preproduct;
            P_len++;
        }
    }
    
    // set primes and exponent arrays for L
    if( init_LofP == 1 ) { L_len = 0; }
    else
    {
        L_len = 1;
        L_distinct_primes[0] = 2;
        L_exponents[0] = __builtin_ctzl(init_LofP );
        init_LofP = ( init_LofP >>  __builtin_ctzl( init_LofP ) );
        // the prime 2 is now taken care of.  Account for odd primes:
        temp = 3;
        while( init_LofP != 1 )
        {
            if( init_LofP % temp == 0 )
            {
                L_exponents[ L_len ] = 0;
                L_distinct_primes[ L_len ] = temp;
                while( init_LofP % temp == 0 )
                {
                    L_exponents[ L_len ]++;
                    init_LofP /= temp;
                }
                L_len++;
            }
            temp += 2;
        }
    }
    len_appended_primes = 0;
}


// assumes prime_stuff is valid and admissible to PP
void Preproduct::appending( Preproduct PP, primes_stuff p )
{
    mpz_mul_ui( P, PP.P, p.prime );
    P_len = PP.P_len + 1;
    std::copy( PP.P_primes,PP.P_primes + PP.P_len, P_primes );
    P_primes[ P_len ] =  p.prime;   
    append_bound = p.prime;
    
    // initialize L to be PP.L and increase only when new factors are seen
    // L = PP.L;
    // if( (PP.L) % (p.prime - 1) == 0 )

    mpz_set( L, PP.L );
    if( mpz_divisible_ui_p( L, ( p.prime - 1 ) )  )
    {
        std::copy( PP.L_distinct_primes,PP.L_distinct_primes + PP.L_len, L_distinct_primes );
        std::copy( PP.L_exponents, PP.L_exponents + PP.L_len, L_exponents );
        L_len = PP.L_len;
    }
    else
    {
        //merge L and p-1
        int i = 0; //counter for PP
        int j = 0; //counter for p
        L_len = 0; //counter for P
        while( i < PP.L_len && j < p.pm1_len )
        {
            // check if they have the same prime
            if( PP.L_distinct_primes[i] == p.pm1_distinct_primes[j] )
            {
                L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
                L_exponents[ L_len ] = std::max( PP.L_exponents[i], p.pm1_exponents[j] );
                //we need to update L if  L_exp < pm1_exp
                if( PP.L_exponents[i] <  p.pm1_exponents[j])
                {
                    for( int L_update = PP.L_exponents[i]; L_update < p.pm1_exponents[j]; L_update++ )
                    {
                        // need to import PP.L_distinct_primes[i] before multplying ?
                        mpz_mul_ui( L, L, PP.L_distinct_primes[i] );
                        //L *= PP.L_distinct_primes[i];
                    }
                }
                i++; j++; L_len++;
            }
            else if( PP.L_distinct_primes[i] < p.pm1_distinct_primes[j] )
            {
                L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
                L_exponents[ L_len ] = PP.L_exponents[i];
                i++; L_len++;
            }
            else
            {
                L_distinct_primes[ L_len ] = p.pm1_distinct_primes[j];
                L_exponents[ L_len ] = p.pm1_exponents[j];
                for( int L_update = 0; L_update < p.pm1_exponents[j]; L_update++ )
                {
                    mpz_mul_ui( L, L, p.pm1_distinct_primes[j] );
                    //L *= p.pm1_distinct_primes[j];
                }
                j++; L_len++;
            }
        }
        while( i < PP.L_len )
        {
            L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
            L_exponents[ L_len ] = PP.L_exponents[i];
            i++; L_len++;
        }
        while( j < p.pm1_len  )
        {
            L_distinct_primes[ L_len ] = p.pm1_distinct_primes[j];
            L_exponents[ L_len ] = p.pm1_exponents[j];
            for( int L_update = 0; L_update < p.pm1_exponents[j]; L_update++ )
            {
                mpz_mul_ui( L, L, p.pm1_distinct_primes[j] );
                // L *= p.pm1_distinct_primes[j];
            }
            j++; L_len++;
        }
    }
    //set appended prime info for further admissibility checks
    len_appended_primes = PP.len_appended_primes + 1;
    // the next inadmissible prime to p.prime is 2*p.prime + 1 or 4*p.prime + 1
    // case is chosen to avoid divisibility by 3
    uint64_t temp_next_inad = 2*p.prime + 1;
    uint16_t temp_mod_3 = 2;
    if( p.prime % 3 == 1 )
    {
        temp_next_inad += 2*p.prime;
        temp_mod_3 = 1; 
    }
    
    if( len_appended_primes == 1 )
    {
        next_inadmissible[0] = temp_next_inad;
        mod_three_status[0] = temp_mod_3;
        appended_primes[0] = p.prime;
    }
    // So, we have an arry in sorted order and we need to insert a new element
    // copy the info from PP into the new product until we find where the new info goes
    // put the new info in
    // resume copying the info from PP if needed
    else
    {
        int i = 0;
        while( PP.next_inadmissible[i] < temp_next_inad && i < PP.len_appended_primes )
        {
            next_inadmissible[i] = PP.next_inadmissible[i];
            mod_three_status[i] = PP.mod_three_status[i];
            appended_primes[i] = PP.appended_primes[i];
            i++;
        }
        next_inadmissible[i] = temp_next_inad;
        mod_three_status[i] = temp_mod_3;
        appended_primes[i] = p.prime; 
        while( i < PP.len_appended_primes )
        {
            next_inadmissible[i+1] = PP.next_inadmissible[i];
            mod_three_status[i+1] = PP.mod_three_status[i];
            appended_primes[i+1] = PP.appended_primes[i];
            i++;
        }
    }
}

// admissibility check with no gcd
// if the while loop is not taken, this will execute with less than 10 instructions
// assumes len_appended_primes > 0
bool Preproduct::is_admissible( uint64_t prime_to_append )
{
  // if this works how we want it to, this while loop will not be entered often
  while( prime_to_append > next_inadmissible[0] )
  {
    // add 2*p or 4*p and avoid divisibility be 3.  Flip the state of mod 3 status.
    next_inadmissible[0] += (appended_primes[0] << mod_three_status[0]);
    mod_three_status[0] = (mod_three_status[0] == 1 ) ? 2 : 1;
    int i= 1;
    // put next_inadmissible in sorted order having increased the first element
    while( i < len_appended_primes && next_inadmissible[ i-1 ] < next_inadmissible[ i ] )
    {
      std::swap( next_inadmissible[ i-1 ], next_inadmissible[ i ] );
      std::swap( appended_primes[ i-1 ], appended_primes[ i ] );
      std::swap( mod_three_status[ i-1 ], mod_three_status[ i ] );
      i++;
    }
  }
  // now prime to append is less than or equal to next_inadmissible
  // return true if it is less than; return false if it is equal
  return ( prime_to_append < next_inadmissible[0] ) ;
}


void Preproduct::CN_search(  )
{
    const uint32_t cache_bound = 150'000;
        
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

    // n will be created in an arithmetic progression with P*L*(possibly lifted small primes)
    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );

    // will need more bases later
    // use bases from the prime divisors of L
    mpz_t base;
    mpz_init( base );

    // storage for the gcd result
    mpz_t gcd_result;
    mpz_init( gcd_result );
    
    mpz_t mpz_prime;
    mpz_init( mpz_prime );
    
    mpz_t small_prime;
    mpz_init( small_prime );
    
    mpz_t lifted_L;
    mpz_init_set( lifted_L, L );
    
    // position i has the truth value of the statement "(2i + 3) is prime"
    std::bitset<256> small_primes{"0010010100010100010000010110100010010100110000110000100010100010010100100100010010110000100100000010110100000010000110100110010010010000110010110100000100000110110000110010010100100110000110010100000010110110100010010100110100110010010110100110010110110111"};
    // fix append_bound to be consistent with the bitset
    uint64_t sieve_index_bound = std::min( (append_bound - 3)/2, (uint64_t) 512);
    
    // remove primes dividing L from the bitset
    uint16_t i = 1;
    while( i < L_len &&  L_distinct_primes[i] < 512 )
    {
        small_primes [ ( L_distinct_primes[i] - 3) /2 ] = 0;
        i++;
    }
    
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    mpz_cdiv_q( cmp_bound, BOUND, PL);
    mpz_t R;
    mpz_init( R );
    
    // store the primes that we append
    std::vector< uint16_t > primes_lifting_L;
    
    // L_lift will store the product of primes and gives the count of spokes on the wheel
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
            mpz_mul_ui( lifted_L, lifted_L, p );
            mpz_cdiv_q( cmp_bound, BOUND, PL);
        }
        prime_index++;
    }

    uint32_t cmp_bound32 = mpz_get_ui( cmp_bound );
    boost::dynamic_bitset<> spoke_sieve( cmp_bound32 );
        
    // This is the wheel and a creates r_star + m*L
    // The wheel is r_star + m*L
    for( uint64_t m = 0; m < L_lift; m ++)
    {
        bool enter_loop = true;
        // checks to make sure that r_star + m*L is not divisible by lifted_primes
        for( auto p : primes_lifting_L )
        {
            enter_loop = ( enter_loop && ( mpz_divisible_ui_p( r_star, p ) == 0 ) );
        }
        
        // This is a spoke.  It sieves on k-values numbers of the form (r_star + m*L) + k*L_lift*L
        if( enter_loop )
        {
            //sieve a spoke
            spoke_sieve.reset();
            uint16_t sieve_prime_index = prime_index;
            while( sieve_prime_index < sieve_index_bound )
            {
                // r_star + k*(L_lift*L) = 0 mod p
                // implies k = -r*( L_lift*L)^{-1} mod p
               if( small_primes[ sieve_prime_index ] )
                {
                    uint32_t p = 2*sieve_prime_index + 3;
                    mpz_set_ui( small_prime, p );
                    // n is being used as a temporary variable in the following 5 lines
                    mpz_invert( n, lifted_L, small_prime );     // n has (L_lift*L)^{-1} mod p
                    mpz_neg( n, n);                             // n has -(L)^{-1} mod p
                    mpz_mul( n, n, r_star );                    // muliply by r_star
                    mpz_mod( n, n, small_prime );               // reduce modulo p
                    uint64_t k = mpz_get_ui( n );

                    while( k < cmp_bound32 )
                    {
                        spoke_sieve[k] = 1;
                        k += p;
                    }
                }
                sieve_prime_index++;
            }
            
            // spoke has been sieved, so only do modular exponentiations on valid places
            // n is initialized correctly here
            mpz_mul( n, P, r_star);
            uint32_t k = 0;
            while( mpz_cmp( n , BOUND ) < 0 )
            {
                if( spoke_sieve[k] == 0 )
                {
                    mpz_powm( fermat_result,  base2,  n, n); // 2^n mod n
                    if( mpz_cmp( fermat_result, base2 ) == 0 )  // check if 2 = 2^n mod n
                    {
                        mpz_powm( fermat_result,  base3,  n, n); // 3^n mod n
                        if( mpz_cmp( fermat_result, base3 ) == 0 )  // check if 3 = 3^n mod n
                        {
                            mpz_divexact( R, n, P);
                            CN_factorization( n, R );
                            // gmp_printf( "n = %Zd = %Zd * %Zd is a base-2 and base-3 Fermat psp. \n", n, P, R);
                        }
                    }
                }
                k++;
                mpz_add( n, n, PL);
            }
        }
        mpz_add( r_star, r_star, L);
    }
       
    mpz_clear( R );
    mpz_clear( small_prime );
    mpz_clear( cmp_bound );
    mpz_clear( r_star );
    mpz_clear( n );
    mpz_clear( PL );
    mpz_clear( base2 );
    mpz_clear( base3 );
    mpz_clear( fermat_result );
    mpz_clear( lifted_L );
        
}

bool Preproduct::appending_is_CN( std::vector< uint64_t >&  primes_to_append )
{
    mpz_t P_temp;
    mpz_t L_temp;
    mpz_t gcd_result;

    mpz_init_set( P_temp, P );
    mpz_init_set( L_temp, L );
    mpz_init_set_ui( gcd_result, 1);
    
    bool return_val = true;

    // not meaningful optimization availalbe
    // break the loop(s) when return_val is detected to be false
    for( auto app_prime : primes_to_append )
    {
        for( int i = 0; i < P_len; i++ )
        {
            return_val = ( return_val && ( app_prime % P_primes[i] != 1 ) );
        }
        // update the preproduct by this prime
        mpz_mul_ui( P_temp, P_temp, app_prime );
        uint64_t temp = app_prime - 1;
        // compute LCM( L, p-1 ) = (L / gcd( L, p-1 ) )*(p-1)
        mpz_gcd_ui( gcd_result, L_temp, temp );
        mpz_divexact( L_temp, L_temp, gcd_result );
        mpz_mul_ui( L_temp, L_temp, temp );
    }

    mpz_sub_ui( P_temp, P_temp, 1);
    return_val = return_val && mpz_divisible_p( P_temp, L_temp );
       
    mpz_clear( P_temp );
    mpz_clear( L_temp );
    mpz_clear( gcd_result );
    
    return return_val;
}

std::vector< primes_stuff > Preproduct::primes_admissible_to_P( )
{
    std::vector< primes_stuff > return_vector;
    
    // there are 5761455 primes less than 10^8
    // we set aside the appropriate space
    // do this or let the space be dynamically allocated based on the computed prime_bound?
    return_vector.reserve( PRIME_COUNT );
    
    // a different way to do the below would be to
    // test if P > 10^8 first
    // only when P > 10^8 would prime_bound have a value less than
    mpz_t sqrt_P;
    mpz_init( sqrt_P );
    mpz_sqrt( sqrt_P, P );
    
    int64_t prime_bound;
    mpz_export( &prime_bound, 0, 1, sizeof(uint64_t), 0, 0, sqrt_P);
    prime_bound = SQRT_BOUND / prime_bound;
    prime_bound = std::min( prime_bound, (int64_t) DEFAULT_MAX_PRIME_BOUND );
    
    /*
     initialize some kind of factor sieve up to prime_bound
     use primes dividing P for admissibility checks
     e.g. for q in the factor sieve make sure 1 != q mod p for each p dividing P
     write each admissible prime into the prime_stuff struct formate
     push_back said prime_stuff
     */
    
    mpz_clear( sqrt_P );
    return return_vector;
}

// Check that lambda(P) divides (P-1)
// consider changing this to int-type return matching how gmp returns
// and have the same return standard as gmp
bool Preproduct::is_CN( )
{
    bool return_val;
    mpz_sub_ui( P, P, 1);
    return_val = mpz_divisible_p( P, L );
    mpz_add_ui( P, P, 1);
    return return_val;
}

/* Factor n, storing the unique prime factors in the associated parameter.
Assume that we already know that n is a base-2 and base-3 Fermat pseudoprime.
Use the strong test to those two bases to split.  Return true if fully factored, false if not
*/

/* Factor n, storing the unique prime factors in the associated parameter.
   Technique is the Fermat method: if n passes the fermat test to a base b, the associated 
   strong test can be used to split.  If n fails a fermat test, quit and return false to signify not carmichael.
*/
bool Preproduct::fermat_factor(uint64_t n, std::vector<uint64_t>& prime_factors)
{
    // bases checked will start at 2, then proceed through odd numbers.
    mpz_t base;
    mpz_init( base );
    mpz_set_ui( base, 2 );

    // storage for the gcd result
    mpz_t gcd_result1;
    mpz_init( gcd_result1 );
    mpz_t gcd_result2;
    mpz_init( gcd_result2 );

    // storage for the strong result, along with its +/- 1 algebraic pieces
    mpz_t strong_result;
    mpz_init( strong_result );
    mpz_t strong_plusone;
    mpz_init( strong_plusone );
    mpz_t strong_minusone;
    mpz_init( strong_minusone );

    // queue for composite factors to be split further.  Starts with n in it
    // and vars for the factors pulled off the queue
    std::queue<uint64_t> composite_factors;
    composite_factors.push(n);
    
    mpz_t r_factor;
    mpz_init( r_factor );
    uint64_t temp;

    bool is_fermat_psp;

    // loop until no more composite factors, or until a fermat test fails
    do
    {
        //testing
        //gmp_printf ("base = %Zd \n", base);
        
        // Fermat test, with b^( (n-1) / 2^e ) stored in strong_result
        is_fermat_psp = fermat_test( n, base, strong_result );
        
        // this conditional is not expected to be entered
        // so the do-while loop is not expected to be invoked
        // most numbers are not Fermat pseudoprimes
        if( is_fermat_psp )
        {
          int start_size = composite_factors.size();
          // use a for loop to go through all factors that are currently in the queue
          for( int j = 0; j < start_size; j++ )
          {
            // get element out of queue and put into mpz_t
            // first time through, this is just r_factor will have the value of r_star
            temp = composite_factors.front();
            composite_factors.pop();
            
            mpz_set_ui( r_factor, temp );

            // loop over all algebraic factors: b^((n-1)/2^k) +/- 1 for e <= k <= 0
            while( mpz_cmp( strong_result, base ) != 0 )
            {
                // check gcd before prime testing
                // strong_plusone holds b^((n-1)/2^k) + 1, strong_minusone holds b^((n-1)/2^k) - 1
                mpz_add_ui( strong_plusone, strong_result, 1 );
                mpz_sub_ui( strong_minusone, strong_result, 1 );
                mpz_gcd( gcd_result1, strong_plusone, r_factor );
                mpz_gcd( gcd_result2, strong_minusone, r_factor );
    
                //gmp_printf ("GCD %Zd results from %Zd and %Zd \n", gcd_result, strong_result, r_factor);
                  
                // check that gcd_result has a nontrivial divisor of n
                // could probably be a check on result1 = +/- 1 mod n
                // before computing the gcd
                if( mpz_cmp(gcd_result1, r_factor) < 0 && mpz_cmp_ui(gcd_result1, 1) > 0 )
                {
                  // will need to add a check about a lower bound on these divisors
                  mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result1);      
                  ( mpz_probab_prime_p( gcd_result1, 0 ) == 0 ) ? composite_factors.push( temp ) : prime_factors.push_back( temp );
                  mpz_divexact( gcd_result1, r_factor, gcd_result1 );
                  mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result1);
                  ( mpz_probab_prime_p( gcd_result1, 0 ) == 0 ) ? composite_factors.push( temp ) : prime_factors.push_back( temp );
                }

                // check gcd( b^((n-1)/2^k) - 1, r_factor ) to see if a split occurred
                else if( mpz_cmp(gcd_result2, r_factor) < 0 && mpz_cmp_ui(gcd_result2, 1) > 0 )
                {
                  // will need to add a check about a lower bound on these divisors
                  mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result2);      
                  ( mpz_probab_prime_p( gcd_result2, 0 ) == 0 ) ? composite_factors.push( temp ) : prime_factors.push_back( temp );
                  mpz_divexact( gcd_result2, r_factor, gcd_result2 );
                  mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result2);
                  ( mpz_probab_prime_p( gcd_result2, 0 ) == 0 ) ? composite_factors.push( temp ) : prime_factors.push_back( temp );
                }
                else // r_factor was not factored, so it is prime or composite
                {
                  ( mpz_probab_prime_p( r_factor, 0 ) == 0 ) ? composite_factors.push( temp ) : prime_factors.push_back( temp );
                }

                // square strong_result, this turns b^((n-1)/2^k) into b^((n-1)/2^(k+1))
                mpz_mul( strong_result, strong_result, strong_result );
            }
          }
          // check if Carmichael, or put that somewhere else?
        }
        
        // get next Fermat base.  If 2, add 1.  If not 2, add 2 to get next odd
        if( mpz_cmp_ui(base, 2) == 0) mpz_add_ui( base, base, 1);
        else mpz_add_ui( base, base, 2 );
        
        // do it again if
        // the number is a Fermat psp and
        // R_composite queue is not empty
        // if i == L_len, we should probably output or factor directly 
            // could be some strange multi-base Fermat pseudoprime - very rare?
        // room for improvement here
    }
    while( is_fermat_psp && !composite_factors.empty() );

    // need to deallocate the mpz_t vars
    mpz_clear( base );
    mpz_clear( gcd_result1 );
    mpz_clear( gcd_result2 );
    mpz_clear( strong_result );
    mpz_clear( r_factor );
    mpz_clear( strong_plusone );
    mpz_clear( strong_minusone );
    
    // if while loop ended because not a fermat psp, result false
    if( !is_fermat_psp ) return false;
    else return true;

}

/* Check whether n is a Fermat pseudoprime to the base b (via b^n = b mod n).  Returns bool with this result.
   Additionally, sets strong_result variable to b^((n-1)/2^e)
   Notes this function returns true for prime n.
*/
bool Preproduct::fermat_test(uint64_t& n, mpz_t& b, mpz_t& strong_result)
{
    // convert n to mpz
    mpz_t n_as_mpz;
    mpz_init( n_as_mpz );
    mpz_set_ui( n_as_mpz, n );
    
    // create a variable for n-1, then compute the largest power of 2 that divides n-1
    mpz_t nminus;
    mpz_init( nminus );
    mpz_sub_ui( nminus, n_as_mpz, 1);
    // counts the number of 0's that terminate in nminus, i.e. e such that 2^e || n-1
    uint32_t exp_on_2 = (uint16_t) mpz_scan1( nminus , 0);
    uint32_t pow_of_2 = ( 1 << exp_on_2 );
    // QUESTION: should pow_of_2 be larger than 32 bits?
    // ANSWER:  no.  Indeed, a char would work.  everything in this computation can fit in an uint128_t
    // which means that exponent is bounded by 128 which is a 7 bit number.
    
    
    // stores the exponent we will apply to the base
    mpz_t strong_exp;
    mpz_init( strong_exp );

    // stores the result b^(n-1) mod n, i.e. the Fermat exponent
    mpz_t fermat_result;
    mpz_init( fermat_result );
    
    // set up strong base:  truncated divsion by 2^e means the exponent holds (n-1)/(2^e)
    mpz_tdiv_q_2exp( strong_exp, nminus, exp_on_2 );
    // we use prime divisors of L as the Fermat bases
    mpz_powm( strong_result,  b,  strong_exp, n_as_mpz); // b^( (n-1)/(2^e) ) mod n
    mpz_powm_ui( fermat_result,  strong_result, pow_of_2, n_as_mpz); // b^( (n-1)/(2^e)) )^(2^e) = b^(n-1)
    mpz_mul( fermat_result, fermat_result, b );  // b^n
    
    //std::cout << "strong_exp = " << mpz_get_str(NULL, 10, strong_exp) << " strong_result = " << mpz_get_str(NULL, 10, strong_result);
    //std::cout << " fermat_result = " << mpz_get_str(NULL, 10, fermat_result) << "\n";

    bool is_psp = mpz_cmp( fermat_result, b ) == 0;

    // deallocating mpz vars created in this function
    mpz_clear( n_as_mpz );
    mpz_clear( nminus );
    mpz_clear( strong_exp );
    mpz_clear( fermat_result );
    
    return is_psp;   
}


// compare with:
// https://github.com/ashallue/tabulate_car/blob/master/LargePreproduct.cpp#L439C1-L500C2
// see section 5.3 of ANTS 2024 work
// assume that B/PL is "big"
// if the prime to be found is of the form r_star + k*L for some relatively small value of k
// call the CN_search method instead (might have to filter out composite completions)
// we do this "unbounded" - if P^2 L > B, it could produce a CN exceeding B
// in which case, it should be discarded after the fact or checked in the inner-loop
void Preproduct::completing_with_exactly_one_prime()
{
    // set up scaled problem:
    mpz_t R;
    mpz_init( R );
    
    mpz_t script_R;
    mpz_init( script_R );
    mpz_invert( script_R, P, L );
    mpz_sub_ui( script_R, script_R, 1 );
    
    //scaled problem in terms of the gcd of r_star - 1 and L
    mpz_t g;
    mpz_init( g );
    mpz_gcd( g, script_R, L );
    
    // the next line now has script_R holding R1 as described in Section 5.3.1 of ANTS 2024 paper
    mpz_divexact( script_R, script_R, g );
    
    // script_P = (P-1)/g
    mpz_t script_P;
    mpz_init_set( script_P, P);
    mpz_sub_ui( script_P, script_P, 1);
    mpz_divexact( script_P, script_P, g);
    
    // script_L = L/g
    mpz_t script_L;
    mpz_init_set( script_L, L);
    mpz_divexact( script_L, script_L, g);
    
    mpz_t div_bound;
    mpz_init_set( div_bound, script_P );
    mpz_sqrt( div_bound, div_bound );
    
    mpz_t divisor;
    mpz_set( divisor, script_R );
    
    /*
     while( script_R + k*script_L < sqrt( script_P ) )
     {
     
     test if g*(script_R + k*script_L) + 1 = r_star + k*L is prime
     if it is, test if n = P * (r_star + k*L ) is a CN
     update.  In terms of script_R or r_star:
     script_R += script_L
     r_star += L
     if the bounds are done as above, no need to explicitly keep track of k
     the bounds are loop-invariant, so we can compute the maximal value of k outside loop and explicitly track k if we want to
     }
     
     
     */
    
    while( mpz_cmp( divisor, div_bound) <= 0 )
    {
        mpz_mul( R, divisor, g);
        mpz_add_ui( R, R, 1);
        if( mpz_probab_prime_p( R, 0 ) != 0 )
        {
            // we might do a bounds check (?)
            // test that P*R is a CN
        }
        mpz_add( divisor, divisor, script_L );
    }
    
    
    
    // set R to be R_2 in section 5.3.2
    mpz_invert( script_R, script_R, script_L);
    mpz_mul( script_R, script_P, script_R );
    mpz_mod( divisor, script_R, script_L );
    
    
    /*
     while(  R_2 + k*script_L < sqrt( script P ) )
     {
     test if R_2 + k*script_L exactly divides (P-1)
     if so, define R = (P-1)/(R_2 + k*script_L) + 1, check if R is prime, and if so, check if PR is a CN with Korselt
     update the arithmetic progression.
     script_R += script_L
     k++
     }
     */
    
    while( mpz_cmp( divisor, div_bound) <= 0 )
    {
        if( mpz_divisible_p( script_P, divisor ) )
        {
            mpz_divexact( R, script_P, divisor);
            mpz_mul( R, R, g);
            mpz_add_ui( R, R, 1);
            if( mpz_probab_prime_p( R, 0 ) != 0 )
            {
                // we might do a bounds check (?)
                // test that P*R is a CN
            }
        }
        mpz_add( divisor, divisor, script_L );
    }
    
    mpz_clear( divisor );
    mpz_clear( div_bound );
    mpz_clear( script_R );
    mpz_clear( script_L );
    mpz_clear( script_P );
    mpz_clear( g );
    mpz_clear( R );

}


// while R is passed as mpz_t type
// the assumption is that R < 2^64
// NEED TO DO - incorporate append bound, see comments below
bool Preproduct::CN_factorization(mpz_t& n, mpz_t& R)
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
        // could compare the least element in the sort
        std::sort ( R_prime_factors.begin(), R_prime_factors.end() );
        
        if( appending_is_CN( R_prime_factors ) )
        {
            std::cout<< "          THIS IS A CARMICHAEL NUMBER     " << std::endl;
            gmp_printf ("%Zd = ", n );
            for( int i = 0; i < P_len; i++ )
                std::cout << " " << P_primes[ i ] ;
            for( int i = 0; i < R_prime_factors.size(); i ++ )
                std::cout << " " << R_prime_factors[i];
            std::cout << std::endl;
        }
    }
        
    mpz_clear( fermat_result );
    mpz_clear( gcd_result );
    mpz_clear( base );
    mpz_clear( odd_part );
    
    return is_fermat_psp;
}



int main(void) {
    
    Preproduct P0;
    // P0.initializing( 599266767, 890750, 991 );
    P0.initializing( 6682828353, 2289560, 13 );
    
    std::cout << "LP for this is " << sizeof( unsigned long int ) << std::endl;
    
    std::cout << "Initializing P : " ;
    gmp_printf ("%Zd = ", P0.P );
    
    for( int i = 0; i < ( P0.P_len - 1 ); i++)
    {
        std::cout << P0.P_primes[i] << " * "  ;       
    }
    std::cout << P0.P_primes[P0.P_len - 1 ] << std::endl ;  
    
    std::cout << "Initializing Lambda : " ;
    gmp_printf ("%Zd = ", P0.L );
    
    for( int i = 0; i < ( P0.L_len - 1 ); i++)
    {
        std::cout << P0.L_distinct_primes[i] << " ^ "  << P0.L_exponents[ i ] << " * "  ;       
    }
    std::cout << P0.L_distinct_primes[P0.L_len - 1] << " ^ "  << P0.L_exponents[ P0.L_len - 1 ] << std::endl ;  

   
    /*  A good test case
    Preproduct P1;
    P1.initializing( 515410417841, 115920, 50 );
    P1.CN_search();
    */
     
    // the triple 132582235, 91872, 941 is in the working jobs with the new rule
    // this means we have between 2 and 6 primes to append to find a CN
    // for all primes in ( 941, (B/P)^(1/3) ) = (941, 196112) we append a prime and do CN_search
    // for all primes in ( 196112, (B/P)^(1/2) ) = ( 196112, 86847502) we append a prime and look for exactly one prime
    // this represents the first two CN_xearches invoked by the above:
    
    
    uint64_t P = 132582235;
    uint64_t L = 91872;
    
    uint64_t P1;
    uint64_t L1;
    
    uint64_t p64 = 947;
    mpz_t p ;
    mpz_init_set_ui( p, p64 );
    
   
    // this example is *not* ideal, we refactor P by trial division every time
    // the dumb_appending isn't working
    // I dont' know why - time to debug
    Preproduct P_testing;
    
    while( p64 < 196112 )
    {
        mpz_nextprime( p, p);
        p64 = mpz_get_ui( p );
        if( 1 != p64 % 5 && 1 != p64 % 17 && 1 != p64 % 23 && 1 != p64 % 73 && 1 != p64 % 929 )
        {
           // std::cout << "I found an admissible prime " << p64 << std::endl;
            P1 = P*p64;
            L1 = L*((p64-1)/std::gcd( L , p64-1 ) );
            P_testing.initializing ( P1, L1, p64 );
            P_testing.CN_search();
            
        }
    }
    /*
     This run completes in 4 minutes and 30 seconds on thomas.butler.edu and finds these:
     433493717815335774335905 =  5 17 23 73 929 2377 7129 10267 18793
     72425097332690148535105 =  5 17 23 73 929 5347 577 673 263089
     433493717815335774335905 =  5 17 23 73 929 7129 2377 10267 18793
     433493717815335774335905 =  5 17 23 73 929 10267 2377 7129 18793
     433493717815335774335905 =  5 17 23 73 929 18793 2377 7129 10267
     725906640592907462305 =  5 17 23 73 929 21577 643 394633
     
     The fact that these are found duplicated only indicates our failure to implement the append_bound-related checks in CN_factorization
     */
   
    
    
    mpz_clear( p );
    
    // P0.CN_search(1873371784);
    //P0.CN_search(149637241475922);
    
    return 0;
}

